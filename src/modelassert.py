#!/usr/bin/env python

""" metano Model Assertion Monitor

This script reads in a metabolic model from a reaction file and a scenario file,
performs FBA, FVA, and split-ratio analysis and checks the results against a
list of assertions. All assertions that do not hold are reported.

It must never be used in a non-safe environment (such as in a server-based
application), because eval() is used to evaluate the assertions.

Allowed symbols are as follows:
<reaction_name>      -- flux through reaction in FBA solution
min(<reaction_name>) -- FVA minimum of flux through reaction
max(<reaction_name>) -- FVA maximum        - '' -
<metabolite_name>    -- total flux going through metabolite (computed as influx)
inratio(<metabolite_name>, <reaction_name>)
                     -- fraction of metabolite produced via reaction
outratio(<metabolite_name>, <reaction_name>)
                     -- fraction of metabolite consumed via reaction


This file is part of metano.
Copyright (C) 2010-2019 Alexander Riemer, Julia Helmecke
Braunschweig University of Technology,
Dept. of Bioinformatics and Biochemistry

metano is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

metano is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with metano.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import range
from builtins import object
from metano.defines import FbaParam, COPYRIGHT_VERSION_STRING
from metano.fba import OptionParser, FbAnalyzer
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.fva import FvAnalyzer
import os
import re


class ModelWatcher(object):
    """ Class for plausibility analysis of metabolic models
    """

    commentSign = '#'

    def __init__(self, model=None, fvaTolerance=.95):
        """ initialize watcher class

        Keyword arguments:

        model        -- MetabolicModel to be checked
        fvaTolerance -- objective function threshold for FVA
        """
        self.setModel(model)
        self.clearAssertions()
        self.fvaTolerance = fvaTolerance

    def clearAssertions(self):
        self.assertions, self.assertions_orig = [], []

    def setModel(self, model):
        """ set the MetabolicModel and reset flux and related variables
        """
        self.model = model
        self.flux = MetabolicFlux()
        self.metFlux, self.inRatios, self.outRatios = {}, {}, {}
        self.minFlux, self.maxFlux = {}, {}

    def readAssertionsFromFile(self, filename):
        """ parse the assertion file identified by filename
            -- wrapper for readAssertionsFromFileHandle()
        """
        with open(filename) as f:
            return self.readAssertionsFromFileHandle(f)

    def readAssertionsFromFileHandle(self, f):
        """ read the assertion file given as a file object

        Modified member variables:

        assertions_orig -- list of assertions in original format
        assertions      -- list of Python statements to be evaluated

        Comments and empty (or all-space) lines are removed
        """
        s = f.read()

        # Remove comments - everything from comment sign ('#') to end of line
        commentPos = 0
        while commentPos >= 0:
            commentPos = s.find(self.commentSign, commentPos)
            if commentPos >= 0:
                eolPos = s.find('\n', commentPos)
                if eolPos < 0:
                    # Last line - delete everything starting from comment sign
                    s = s[:commentPos]
                    commentPos = -1
                else:
                    # Excise comment
                    s = s[:commentPos]+s[eolPos:]

        if s.find(':') >= 0:
            raise SyntaxError('Illegal character found (:)')

        self.assertions_orig = s.split('\n')

        # Replace reaction and metabolite names with unique non-matching IDs
        # including a type specifier (r or m)

        idDict = dict(list(zip(list(self.model.reactionDict.keys()), "r"*len(self.model))) +
                      list(zip(list(self.model.metaboliteDict.keys()),
                               "m"*len(self.model.metaboliteDict))))

        # Find non-overlapping longest matches (greedy) of reaction and
        # metabolite names and replace with unique non-matching IDs
        counter = 0
        replaceDict, replaceDictRev = {}, {}
        for ID in sorted(list(idDict.keys()), key=len, reverse=True):
            pos = 0
            while pos >= 0:
                pos = s.find(ID, pos)
                if pos >= 0:
                    if ID in replaceDict:
                        replaceId = replaceDict[ID]
                    else:
                        # Generate new ID
                        replaceId = ":%u%c:" % (counter, idDict[ID])
                        replaceDict[ID] = replaceId
                        replaceDictRev[replaceId] = ID
                        counter += 1

                    s = s[:pos]+replaceId+s[pos+len(ID):]
                    pos = pos+len(replaceId)+1

        reaString = r"(:\d+r:)"
        metString = r"(:\d+m:)"
        inratioPattern = re.compile(r"inratio\s*\(\s*"+metString+r"\s*,\s*" +
                                    reaString+r"\s*\)")
        outratioPattern = re.compile(r"outratio\s*\(\s*"+metString+r"\s*,\s*" +
                                     reaString+r"\s*\)")
        minPattern = re.compile(r"min\s*\(\s*"+reaString+r"\s*\)")
        maxPattern = re.compile(r"max\s*\(\s*"+reaString+r"\s*\)")
        reaPattern = re.compile(reaString)
        metPattern = re.compile(metString)

        # Replace inratio/outratio expressions with references to self.inRatios/
        #                                                         self.outRatios
        s = inratioPattern.sub(lambda x: "self.inRatios[%r].get(%r, 0.)" %
                               (replaceDictRev[x.group(1)],
                                replaceDictRev[x.group(2)]), s)
        s = outratioPattern.sub(lambda x: "self.outRatios[%r].get(%r, 0.)" %
                                (replaceDictRev[x.group(1)],
                                 replaceDictRev[x.group(2)]), s)

        # Replace min/max expressions with references to self.minFlux/
        #                                                self.maxFlux
        s = minPattern.sub(lambda x: "self.minFlux[%r]" %
                           replaceDictRev[x.group(1)], s)
        s = maxPattern.sub(lambda x: "self.maxFlux[%r]" %
                           replaceDictRev[x.group(1)], s)

        # Replace reaction identifiers with references to self.flux
        s = reaPattern.sub(lambda x: "self.flux[%r]" %
                           replaceDictRev[x.group(0)], s)

        # Replace metabolite identifiers with references to self.metFlux
        s = metPattern.sub(lambda x: "self.metFlux[%r]" %
                           replaceDictRev[x.group(0)], s)

        self.assertions = s.split('\n')

        # Remove blank lines (cannot be evaluated)
        for i in range(len(self.assertions)-1, -1, -1):
            if self.assertions[i] == "" or self.assertions[i].isspace():
                del self.assertions[i]
        for i in range(len(self.assertions_orig)-1, -1, -1):
            if (self.assertions_orig[i] == "" or
                    self.assertions_orig[i].isspace()):
                del self.assertions_orig[i]

    def checkSyntax(self):
        """ check assertions for syntax errors without actually performing any
            potentially time-consuming analyses

        This function prints a message to the console for every assertion that
        contains a syntax error.

        Returns: True if all assertions can be evaluated, False if error occurs
        """
        # Backup member variables and re-initialize with dummy values to make
        # sure that all references exist
        bakflux, bakMetFlux = self.flux, self.metFlux
        bakInRatios, bakOutRatios = self.inRatios, self.outRatios
        bakMinFlux, bakMaxFlux = self.minFlux, self.maxFlux

        self.flux = MetabolicFlux(self.model, [0.]*len(self.model))
        self.metFlux = dict(list(zip(self.model.metabolites,
                                     [0.]*len(self.model.metabolites))))
        self.inRatios = dict(list(zip(self.model.metabolites,
                                      [{}]*len(self.model.metabolites))))
        self.outRatios = dict(list(zip(self.model.metabolites,
                                       [{}]*len(self.model.metabolites))))
        self.minFlux = dict(list(zip(list(self.model.reactionDict.keys()),
                                     [0.]*len(self.model))))
        self.maxFlux = dict(list(zip(list(self.model.reactionDict.keys()),
                                     [0.]*len(self.model))))
        excepted = []
        for i in range(len(self.assertions)):
            try:
                eval(self.assertions[i])
            except (SyntaxError, NameError):
                excepted.append(repr(self.assertions_orig[i]))

        # Restore member variables to original values
        self.flux, self.metFlux = bakflux, bakMetFlux
        self.inRatios, self.outRatios = bakInRatios, bakOutRatios
        self.minFlux, self.maxFlux = bakMinFlux, bakMaxFlux

        if excepted:
            print ("The following assertions cannot be evaluated due to errors"
                   ":\n  " + "\n  ".join(excepted))
#  The following was removed because the Python error messages are misleading.
#            for a in excepted:
#                print "\n  %r\n  %s: %s\n" % (a[0], a[1], a[2])
            return False
        return True

    def run(self, fbaParams, syntaxOk=True, printExcepted=False):
        """ check assertions against results of FBA, FVA, & split-ratio analysis

        This function performs FBA, split-ratio analysis, and FVA and evaluates
        the assertions. It enumerates all assertions that do not hold.

        Keyword arguments:
        fbaParams     -- objective function and solver for FBA
        syntaxOk      -- if False, previously performed syntax check has failed
        printExcepted -- if False, syntax errors in assertions are not reported
        """

        # 1. Perform flux balance analysis on the metabolic network

        fba = FbAnalyzer(fbaParams.solver)
        self.flux = fba.runOnModel(self.model, fbaParams)[1]
        if len(self.flux) == 0:
            print ("No FBA solution was obtained. The optimization problem is "
                   "either unbounded or infeasible.")
            return

        # 2. Compute metabolite fluxes and split ratios

        splitRatios = self.flux.computeAllSplitRatios(self.model)
        self.metFlux, self.outRatios, self.inRatios = {}, {}, {}
        for met in splitRatios:
            outRatios, inRatios = splitRatios[met]
            self.metFlux[met] = sum(inRatios[rea][1] for rea in inRatios)
            self.outRatios[met] = dict((rea, outRatios[rea][0])
                                       for rea in outRatios)
            self.inRatios[met] = dict((rea, inRatios[rea][0])
                                      for rea in inRatios)

        # 3. Perform flux variability analysis on the model

        fva = FvAnalyzer("default")  # always use GLPK (=> fastFVA)
        minmax = fva.runOnModel(self.model, fbaParams, self.fvaTolerance,
                                self.flux)[0]

        self.minFlux, self.maxFlux = {}, {}
        for i in range(len(self.model)):
            self.minFlux[self.model.reactions[i].name], self.maxFlux[
                self.model.reactions[i].name] = minmax[i]

        # 4. Check assertions

        failed = []
        excepted = []
        for i in range(len(self.assertions)):
            try:
                if not eval(self.assertions[i]):
                    failed.append(repr(self.assertions_orig[i]))
            except (SyntaxError, NameError):
                excepted.append(repr(self.assertions_orig[i]))

        if printExcepted and excepted:
            print ("The following assertions could not be evaluated due to "
                   "errors:\n  " + "\n  ".join(excepted))
        if failed:
            print("The following assertions failed:\n  " + "\n  ".join(failed))
        else:
            if not syntaxOk or (printExcepted and excepted):
                print("All other assertions hold.")
            else:
                print("All assertions hold.")


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "Model assertion watcher\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile", help="read "
                      "metabolic model from the given reaction FILE",
                      metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile", help="use the "
                      "given scenario FILE for model analysis ", metavar="FILE")
    parser.add_option("-a", "--assertions", dest="assertionFile", help="read "
                      "assertions from FILE",  metavar="FILE")
    parser.add_option("-t", "--tolerance", dest="tolerance", type="float",
                      help="FVA: tolerance for objective function (0 < VALUE < "
                      "1; default: .95)", metavar="VALUE")
    parser.set_defaults(tolerance=.95)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-a")

    if options.tolerance <= 0. or options.tolerance > 1.:
        print("Error: Tolerance must be in interval (0, 1]")
        exit()

    # 2. Create MetabolicModel from reaction file

    rparser = ReactionParser()
    model = MetabolicModel()
    try:
        model.addReactionsFromFile(options.reactionFile, rparser)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.reactionFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in reaction file %s:" %
               os.path.basename(options.reactionFile))
        print(strerror)
        exit()

    # 3. Parse scenario file

    model_messages = []
    pparser = ParamParser()
    try:
        # Parse file, get maxmin, name of objective function, and solver name
        maxmin, objStr, solver, numIter, lb, ub = \
            pparser.parse(options.paramFile)
        if solver == "":
            solver = "default"
        fbaParams = FbaParam(solver, maxmin, objStr, numIter)
        # Set flux bounds in model
        model.setFiniteBounds(lb, ub, True, model_messages)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.paramFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in scenario file %s:" %
               os.path.basename(options.paramFile))
        print(strerror)
        exit()
    except ValueError as strerror:
        print(strerror)
        exit()

    # 4. Create ModelWatcher object and parse assertion file

    maw = ModelWatcher(model, options.tolerance)
    try:
        maw.readAssertionsFromFile(options.assertionFile)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.assertionFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in assertion file %s:" %
               os.path.basename(options.assertionFile))
        print(strerror)
        exit()

    # 5. First perform syntax check, then evaluate assertions

    maw.run(fbaParams, maw.checkSyntax())


if __name__ == "__main__":
    main()
