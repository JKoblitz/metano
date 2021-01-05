#!/usr/bin/env python

""" This script checks whether the given flux distribution (FBA solution file)
    conforms to the constraints imposed by the given reaction and/or scenario
    file. It lists all constraints violations. There is zero tolerance, e.g. a
    value of -1E-25 violates an "LB 0." constraint.

Usage:

    to_check_constraints.py <solution-file> [-r reaction-file] [-p param-file]


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
from __future__ import division

from builtins import range
from past.utils import old_div
from metano.fba import OptionParser
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.defines import COPYRIGHT_VERSION_STRING
import os
import sys
from numpy import dot


def checkIrrev(solution, model):
    """ check the signs of the fluxes in 'solution' associated with irreversible
        reactions in the given model

    Returns negViolations, posViolations - {reaction : flux} dicts
    """
    negViolations, posViolations = {}, {}
    for rea in solution:
        if rea not in model.reactionDict:
            continue
        rea_index = model.reactionDict[rea]
        flux = solution[rea]

        if model.reactions[rea_index].lb >= 0. and flux < 0.:
            negViolations[rea] = flux
        elif model.reactions[rea_index].ub <= 0. and flux > 0.:
            posViolations[rea] = flux

    return negViolations, posViolations


def checkBoundsExplicit(solution, lb, ub):
    """ check if flux distribution 'solution' violates any of the given lower
        and upper flux bounds

    Keyword arguments:

    solution -- flux distribution (MetabolicFlux object)
    lb, ub   -- dicts {reaction : bound} of lower and upper bounds

    Returns:  lbViolations, ubViolations - {reaction : flux, bound} dicts
    """
    lbViolations, ubViolations = {}, {}
    for rea in solution:
        flux = solution[rea]

        if rea in lb:
            if flux < lb[rea]:
                lbViolations[rea] = flux, lb[rea]
        if rea in ub:
            if flux > ub[rea]:
                ubViolations[rea] = flux, ub[rea]

    return lbViolations, ubViolations


def checkBoundsByModel(solution, model):
    """ check if the flux distribution 'solution' violates any flux bounds
        defined for the reactions in the given model

    Returns:  lbViolations, ubViolations - {reaction : flux, bound} dicts
    """
    lbViolations, ubViolations = {}, {}
    for rea in solution:
        if rea not in model.reactionDict:
            continue
        rea_index = model.reactionDict[rea]
        flux = solution[rea]
        lb, ub = model.reactions[rea_index].lb, model.reactions[rea_index].ub

        if flux < lb:
            lbViolations[rea] = flux, lb
        elif flux > ub:
            ubViolations[rea] = flux, ub

    return lbViolations, ubViolations


def checkBoundsInSolution(solution):
    """ check if flux distribution 'solution' violates any of the flux bounds
        associated with its reactions

    Returns:  lbViolations, ubViolations - {reaction : flux, bound} dicts
    """
    lbViolations, ubViolations = {}, {}
    for rea in solution:
        flux = solution[rea]
        if rea not in solution.boundsDict:
            continue
        lb, ub = solution.boundsDict[rea]

        if flux < lb:
            lbViolations[rea] = flux, lb
        elif flux > ub:
            ubViolations[rea] = flux, ub

    return lbViolations, ubViolations


def checkMassBalance(solution, model):
    """ check if flux distribution 'solution' violates the mass balance
        constraints defined by the model's stoichiometric matrix

    Returns:
    avgErr, maxPosErr, minNegErr, maxRowMet, minRowMet

    avgErr    -- average deviation from zero over all rows
    maxPosErr -- highest positive deviation from zero
    minNegErr -- lowest negative deviation from zero
    maxRowMet -- row (metabolite name) where maxPosErr occurred
    minRowMet -- row (metabolite name) where minNegErr occurred
    """
    solutionVec = solution.getVecOrderedByModel(model)
    sv = dot(model.getStoichiometricMatrix(), solutionVec)
    svmax = svmin = 0.
    svmax_i = svmin_i = -1
    svsum = 0.
    nMetabolites = len(sv)
    for i in range(nMetabolites):
        svsum += abs(sv[i])
        if sv[i] < svmin:
            svmin = sv[i]
            svmin_i = i
        elif sv[i] > svmax:
            svmax = sv[i]
            svmax_i = i

    return (old_div(svsum, nMetabolites), svmax, svmin, model.metabolites[svmax_i],
            model.metabolites[svmin_i])


def printMassBalanceCheckResult(svavg, svmax, svmin, svmaxMet, svminMet):
    """ print the result of checkMassBalance() - simple call:
        printMassBalanceCheckResult(*checkMassBalance(solution, model))
    """
    print("Mean deviation from zero in S*v: %g" % svavg)
    print ("Maximum pos. deviation from zero in S*v: %g in metabolite %s" %
           (svmax, svmaxMet))
    print ("Maximum neg. deviation from zero in S*v: %g in metabolite %s" %
           (svmin, svminMet))


def checkLinConstraints(solution, model, linConstraints):
    """ check if flux distribution 'solution' violates any of the given linear
        equality and inequality constraints

    Keyword arguments:

    solution       -- flux distribution (MetabolicFlux object)
    model          -- MetabolicModel (for generating flux and coefficient
                                      vectors)
    linConstraints -- list of LinearConstraint objects

    Returns:
    eqErr, ineqErr

    eqErr   -- dict {index in linConstraints : deviation} for each equality
               constraint
    ineqErr -- dict {index in linConstraints : deviation} for each inequality
               constraint - if deviation is negative, constraint is satisfied
    """
    solutionVec = solution.getVecOrderedByModel(model)
    eqs, ineqs = ParamParser.linConstraintsToVectors(linConstraints,
                                                     model.reactionDict)
    eqIndex, ineqIndex = [], []
    for i in range(len(linConstraints)):
        if linConstraints[i].isEq:
            eqIndex.append(i)
        else:
            ineqIndex.append(i)

    eqErr, ineqErr = {}, {}
    if eqs:
        A = [row[0] for row in eqs]
        b = [row[1] for row in eqs]
        dotVec = dot(A, solutionVec)
        for i in range(len(b)):
            eqErr[eqIndex[i]] = abs(dotVec[i]-b[i])

    if ineqs:
        A = [row[0] for row in ineqs]
        b = [row[1] for row in ineqs]
        dotVec = dot(A, solutionVec)
        for i in range(len(b)):
            ineqErr[ineqIndex[i]] = dotVec[i]-b[i]

    return eqErr, ineqErr


def printLinConstraintCheckResult(eqErr, ineqErr, linConstraints):
    if eqErr:
        print("Deviation in additional linear equality constraints:")
        j = 1  # Counter for enumeration of equality constraints
        for i in sorted(eqErr):
            print("%3d: %g (%s)" % (j, eqErr[i], linConstraints[i]))
            j += 1

    if ineqErr:
        conforms = True
        j = 1  # Counter for enumeration of violated inequality constraints
        for i in sorted(ineqErr):
            if ineqErr[i] > 0.:
                if conforms:
                    print("Deviation in linear inequality constraints:")
                    conforms = False
                print("%3d: %g (%s)" % (j, ineqErr[i], linConstraints[i]))
        if conforms:
            print ("The given solution conforms to all linear inequality "
                   "constraints.")


def main():
    # 1. Parse command line

    usage = "Usage: %prog <solution-file> [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile", help="use "
                      "reversibilities from reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use constraints from the given scenario FILE",
                      metavar="FILE")

    options, args = parser.parse_args()

    # 2. Read solution file

    solution = MetabolicFlux()
    try:
        solution.readFromFile(args[0])
    except IndexError:
        print ("Error: No solution file given.\nUsage is\n    " +
               os.path.basename(sys.argv[0]) + " <solution-file> [options]")
        exit()
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(args[0]))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("An error occurred parsing file %s:" %
               os.path.basename(args[0]))
        print(strerror)
        exit()

    if options.reactionFile is None and options.paramFile is None:
        print("Neither reaction nor scenario file given. Nothing to do.")
        exit()

    # 3. Parse reaction file (if given) and check irreversibilities

    if options.reactionFile:
        model = MetabolicModel()
        rparser = ReactionParser()
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

        negViolations, posViolations = checkIrrev(solution, model)
        conforms = not negViolations and not posViolations

        for rea in negViolations:
            print ("Irreversible reaction %s has negative flux (%g)." %
                   (rea, negViolations[rea]))
        for rea in posViolations:
            print ("Flipped irreversible reaction %s has positive flux (%g)." %
                   (rea, posViolations[rea]))
    else:
        conforms = True
        model = None

    # 4. Parse scenario file (if given) and check LB/UB constraints

    if options.paramFile:
        pparser = ParamParser()
        try:
            lb, ub = pparser.parse(options.paramFile)[4:]
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
        linConstraints = pparser.lin_constraints

        lbViolations, ubViolations = checkBoundsExplicit(solution, lb, ub)
        conforms = conforms and not lbViolations and not ubViolations

        for rea in lbViolations:
            print ("%s (%g) violates LB %g." % ((rea,) + lbViolations[rea]))
        for rea in ubViolations:
            print ("%s (%g) violates UB %g." % ((rea,) + ubViolations[rea]))
    else:
        linConstraints = None

    if conforms:
        print ("The given solution does not violate any bounds.")

    if model:

        # 5. Check mass balance constraints

        printMassBalanceCheckResult(*checkMassBalance(solution, model))

        # 6. Check additional equality and inequality constraints

        printLinConstraintCheckResult(*checkLinConstraints(solution, model,
                                                           linConstraints),
                                      linConstraints=linConstraints)


if __name__ == "__main__":
    main()
