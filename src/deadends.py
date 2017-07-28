#!/usr/bin/env python

""" Script for identifying dead ends in metabolic networks

This script identifies dead-end metabolites in a metabolic network, i.e.
metabolites that can either only be produced or only be consumed.

Reactions in which dead-end metabolites occur are assigned a zero flux in
flux-balance analysis, as the net flux through the network must be zero (steady
state condition).

This script wraps function MetabolicModel.findDeadEnds().

With parameter "-a", this script runs a recursive analysis to determine all
reactions that will have zero flux in FBA due to dead ends.
Reactions marked as transporters (identified via reg-ex) are temporarily
activated for the analysis, i.e. the flux bounds for these reactions are
ignored. (For FBA simulations, only a few of these reactions are active at a
time to provide a defined medium.)


This file is part of metano.
Copyright (C) 2010-2017 Alexander Riemer, Julia Helmecke
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

from reactionparser import ReactionParser
from paramparser import ParamParser
from metabolicmodel import MetabolicModel
from fba import OptionParser
from defines import COPYRIGHT_VERSION_STRING
import os, re

# Two floats are considered equal if their difference is less than EPSILON
EPSILON = 1E-8


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "Dead end finder\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="build stoichiometric matrix from reaction FILE",
                      metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use the given scenario FILE", metavar="FILE")
    parser.add_option("-a", "--recursive", action="store_true",
                      dest="recursive", help="recursively predict all fluxes"
                      "that are zero in FBA due to dead ends")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", help="generate extra (debug) output")
    parser.add_option("-t", "--transporter-regex", dest="regex", help="treat "
                      "reactions matching regular expression PATTERN as "
                      "transporters, i.e. temporarily ignore bounds",
                      metavar="PATTERN")
    parser.set_defaults(recursive=False, verbose=False)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")

    # 2. Parse reaction file

    rparser = ReactionParser()
    model = MetabolicModel()
    try:
        model.addReactionsFromFile(options.reactionFile, rparser)
    except IOError, strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.reactionFile))
        print strerror
        exit()
    except SyntaxError, strerror:
        print ("Error in reaction file %s:" %
               os.path.basename(options.reactionFile))
        print strerror
        exit()

    # Get set of transporter reactions matching regular expression
    if options.regex:
        try:
            pattern = re.compile(options.regex)
        except re.error, strerror:
            print ("Error in regular expression: %s. Unable to evaluate." %
                   strerror)
            exit()

        transporters = set(reaction.name for reaction in model if
                           pattern.search(reaction.name))

        if options.verbose:
            print "Ignoring bounds for the following transporters:"
            print ", ".join(sorted(transporters))
    else:
        transporters = set()

    # 3. Parse scenario file

    model_messages = []
    pparser = ParamParser()
    try:
        # Parse file for finite lower and upper bounds
        lb, ub = pparser.parse(options.paramFile)[4:]
        # Set flux bounds in model
        model.setFiniteBounds(lb, ub, True, model_messages)
    except IOError, strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.paramFile))
        print strerror
        exit()
    except SyntaxError, strerror:
        print ("Error in scenario file %s:" %
               os.path.basename(options.paramFile))
        print strerror
        exit()
    except ValueError, strerror:
        print strerror
        exit()

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print '\n'.join(msgs)+'\n'

    # 4. Perform dead-end analysis

    deadEnds, deadReactions = \
        model.findDeadEnds(options.recursive, transporters, options.verbose,
                           EPSILON)

    # 5. Print results

    if deadEnds == []:
        print "All metabolites can be produced or consumed in some reaction."
        exit()

    print ("%u dead-end metabolites (never produced or never consumed):" %
           len(deadEnds))
    for d in sorted(deadEnds):
        print d

    if deadReactions:
        print ("\n%u reactions have zero flux due to dead-end metabolites" %
               len(deadReactions))
        for d in sorted(deadReactions):
            print d


if __name__ == "__main__":
    main()
