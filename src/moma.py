#!/usr/bin/env python

""" metano Minimization of Metabolic Adjustment (MOMA)

This script performs Minimization of Metabolic Adjustment (MOMA) on the
metabolic network (given as a set of reaction equations) and inequality
constraints of a mutant organism against the given FBA solution for the
corresponding wildtype.

It parses a reaction file and a scenario file, formulates the corresponding
linear, quadratic, or other nonlinear optimization problem, and calls an
appropriate solver.


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

from builtins import map
from builtins import range
from builtins import object
from metano.fba import OptionParser, FbAnalyzer
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.momaproblem import MomaProblem
from metano.nonlinearproblem import NonLinearProblem
from metano.startpointiterator import StartPointIterator
from metano.fva import FvAnalyzer
from metano.defines import (SolverStatus, printHistogram, printStatus, FbaParam,
                            COPYRIGHT_VERSION_STRING)
from numpy import array, dot, nan, sqrt, ndarray
from math import exp
import os

_DEFAULT_ALPHA = 1E-6
_DEFAULT_BETA = .15


class MomaAnalyzer(object):
    """ Class for Minimization of Metabolic Adjustment (MOMA)
    """

    def __init__(self, solver="default"):
        self.solver = solver

    @staticmethod
    def getWeightsFromFluxVar(model, fvaMinmax, alpha=_DEFAULT_ALPHA,
                              beta=_DEFAULT_BETA):
        """ compute weights for wMOMA from the given FVA solution

        Keyword arguments:

        model        -- the MetabolicModel
        fvaMinmax    -- dict { reaction : flux minimum, flux maximum }
        alpha, beta  -- parameters for weights computation:
                        w = alpha + exp(-beta*(fvaMax-fvaMin))

        Returns vector of weights (list indexed like model's reactions list)
        """
        weights = []
        for rea in model.getReactionNames():
            mini, maxi = fvaMinmax[rea]
            weights.append(alpha+exp(beta*(mini-maxi)))
        return weights

    def run(self, matrix, lb, ub, wtSolution, eqs=[], ineqs=[], numIter=1,
            weights=None):
        """ construct and solve the quadratic optimization problem

        Keyword arguments:

        matrix     -- the stoichiometric matrix of the metabolic network
        lb         -- list of lower bounds (indexed like matrix columns)
        ub         -- list of upper bounds (indexed like matrix columns)
        wtSolution -- FBA solution for the wildtype
                     (list indexed like matrix columns)
        eqs        -- list of (coefficient vector, right-hand side) pairs for
                      additional equality constraints
        ineqs      -- list of - '' - for additional inequality constraints
        numIter    -- number of iterations of NLP to perform
        weights    -- weight vector for weighted MOMA (None -> perform regular
                      MOMA, else: weight flux i with weights[i])

        Returns:
        distance, solution, status

        distance   -- minimum possible distance from wtSolution with the given
                      matrix & constraints
        solution   -- a flux vector with minimal distance to wtSolution
                     (list indexed like matrix columns)
        status     -- SolverStatus after optimization
        """
        wtVec = array(wtSolution)
        # Construct matrices of original equality and inequality constraints
        Aeq, beq, Aineq, bineq = FbAnalyzer.makeConstraintMatrices(matrix, eqs,
                                                                   ineqs)
        try:
            if MomaProblem.isCvxQpSolver(self.solver):
                ps = MomaProblem(Aeq, beq, Aineq, bineq, lb, ub, self.solver,
                                 wtVec, weights)
            else:
                ps = NonLinearProblem(Aeq, beq, Aineq, bineq, lb, ub,
                                      self.solver)
                if weights is None:
                    ps.setObjective(lambda x: dot(x-wtVec, x-wtVec))
                    ps.setObjGrad(lambda x: 2*(x-wtVec))
                else:
                    if not isinstance(weights, ndarray):
                        weights = array(weights)
                    ps.setObjective(lambda x: dot(x-wtVec, (x-wtVec)*weights))
                    ps.setObjGrad(lambda x: 2*sqrt(weights)*(x-wtVec))
                spIter = StartPointIterator(len(lb), numIter)
                spIter.setRange(-1., 1.)
                ps.setStartPointIterator(spIter)

        except ValueError as strerror:
            print(strerror)
            exit()

        try:
            distance, solution = ps.minimize(gtol=1e-5, maxIter=200,
                                             maxFunEvals=1e5, abstol=1e-6,
                                             feastol=1e-5, reltol=1e-5,
                                             maxiters=100)
        except ValueError as strerror:
            if self.solver == "default":
                print("Default solver reported an error:")
            else:
                print("Solver %s reported an error:" % self.solver)
            print(strerror)
            print("Consider trying a different solver.")
            return nan, [], ps.status
        except KeyError as strerror:
            print(strerror)
            exit()

        if hasattr(ps, "getHistogram"):
            printHistogram(*ps.getHistogram(32))

        return distance, solution, ps.status

    def runOnModel(self, model, wtSolution, linConstraints=[], numIter=1,
                   weights=None, blockedReactions=[]):
        """ construct and solve the quadratic optimization problem
          - this function runs directly on a MetabolicModel and a MetabolicFlux

        Keyword arguments:

        model          -- the MetabolicModel
        wtSolution     -- FBA solution for the wildtype (given as MetabolicFlux)
        linConstraints -- list of LinearConstraint objects
        numIter        -- number of iterations of NLP to perform
        weights        -- weight vector for weighted MOMA (None -> perform
                          regular MOMA, else: weight flux i with weights[i])
        blockedReactions
                       -- remove the given blocked reactions before analysis
                          (faster and gives an optimal solution, as well)

        Returns:
        distance, solution, status, dimReduced

        distance   -- minimum possible distance from wtSolution with the given
                      matrix & constraints
        solution   -- a solution with minimal distance to wtSolution
                     (as MetabolicFlux)
        status     -- SolverStatus after optimization
        dimReduced -- pair (nRows, nColumns) with dimensions of reduced matrix
        """
        if blockedReactions:
            modelRed = model.getSubModelByExcludeList(blockedReactions)
            cbz = model.canBeZero(blockedReactions)
            nonZeroDeadEnds = [blockedReactions[i] for i in
                               range(len(blockedReactions)) if not cbz[i]]
            if nonZeroDeadEnds:
                print ("The following blocked reactions are constrained to a "
                       "non-zero flux:\n  " + "\n  ".join(nonZeroDeadEnds) +
                       "\nThe problem is infeasible.")
                return (nan, MetabolicFlux(), SolverStatus.PRIM_INFEAS,
                        array(modelRed.getStoichiometricMatrix()).shape)

            reactionsRed = set(modelRed.getReactionNames())
            if weights is None:
                weightsRed = None
            else:
                weightsRed = [0.]*len(modelRed)
                for rea in wtSolution:
                    if rea in reactionsRed:
                        weightsRed[modelRed.reactionDict[rea]] = \
                            weights[model.reactionDict[rea]]
                    weightsRed = array(weightsRed)
        else:
            modelRed = model
            weightsRed = weights

        matrix = array(modelRed.getStoichiometricMatrix())
        dimReduced = matrix.shape
        lb, ub = list(map(array, modelRed.getBounds()))
        try:
            eqs, ineqs = ParamParser.linConstraintsToVectors(linConstraints,
                                                             modelRed.reactionDict)
        except ValueError:
            # If any linear constraint is contradictory, return empty solution
            return (nan, MetabolicFlux(), SolverStatus.PRIM_INFEAS,
                    array(modelRed.getStoichiometricMatrix()).shape)

        distance, solution, status = self.run(matrix, lb, ub,
                                              wtSolution.getVecOrderedByModel(
                                                  modelRed),
                                              eqs, ineqs, numIter, weightsRed)
        flux = MetabolicFlux(modelRed, solution)

        # Add removed reactions with flux 0. and original bounds to solution
        if len(modelRed) != len(model) and solution != []:
            reactionsRed = set(modelRed.getReactionNames())
            for rea in model:
                if rea.name not in reactionsRed:
                    flux.fluxDict[rea.name] = 0.
                    flux.boundsDict[rea.name] = (rea.lb, rea.ub)

        return distance, flux, status, dimReduced

    def evalObjFunc(self, wtVec, solution, weights=None):
        if MomaProblem.isCvxQpSolver(self.solver):
            ps = MomaProblem([[0.]*len(wtVec)], lb=[0.]*len(wtVec),
                             ub=[0.]*len(wtVec), solver=self.solver,
                             solution=wtVec, weights=weights)
        else:
            ps = NonLinearProblem([[0.]], lb=[0.], ub=[0.], solver=self.solver)
            if weights is not None:
                if not isinstance(weights, ndarray):
                    weights = array(weights)
                ps.setObjective(lambda x: dot(x-wtVec, (x-wtVec)*weights))
            else:
                ps.setObjective(lambda x: dot(x-wtVec, x-wtVec))
        return ps.eval(solution)


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = ("Minimization of metabolic adjustment\n" +
               COPYRIGHT_VERSION_STRING)
    parser = OptionParser(usage=usage, version=version)
    momaOptions = parser.add_option_group("General MOMA options")
    momaOptions.add_option("-r", "--reactions", dest="reactionFile",
                           help="perform MOMA on the mutant network given by "
                           "the reaction FILE", metavar="FILE")
    momaOptions.add_option("-p", "--parameters", dest="paramFile",
                           help="use the given scenario FILE for the mutant",
                           metavar="FILE")
    momaOptions.add_option("-w", "--wt-solution", dest="wtSolution",
                           help="use the given solution FILE for the wildtype",
                           metavar="FILE")
    momaOptions.add_option("-o", "--output", dest="outputFile",
                           help="write flux distribution computed by MOMA to "
                           "FILE", metavar="FILE")
    momaOptions.add_option("-v", "--reduce-by-fva", dest="redFvaFile",
                           help="use FVA result from FILE to reduce matrix "
                           "(experimental; only works if perturbed solution "
                           "space is a subset of the wildtype solution space)",
                           metavar="FILE")
    wmOptions = parser.add_option_group("Extra options for weighted MOMA")
    wmOptions.add_option("-x", "--wmoma", dest="wMomaFvaFile", help="perform "
                         "weighted MOMA with weights computed from the given "
                         "FVA solution FILE", metavar="FILE")
    wmOptions.add_option("-a", "--alpha", dest="alpha", type="float",
                         help="parameter ALPHA for computation of weights: "
                         "w = alpha + exp(-beta*(fvaMax-fvaMin)) (default: %g)"
                         % _DEFAULT_ALPHA)
    wmOptions.add_option("-b", "--beta", dest="beta", type="float",
                         help="parameter BETA for computation of weights: "
                         "w = alpha + exp(-beta*(fvaMax-fvaMin)) (default: %g)"
                         % _DEFAULT_BETA)
    optOptions = parser.add_option_group("Options for optimization")
    optOptions.add_option("-s", "--solver", dest="solver", help="QP/NLP solver "
                          "to be used for MOMA (default: cvxopt)")
    optOptions.add_option("-i", "--iterations", dest="numIter", type="int",
                          help="number of NLP runs to perform (from different "
                          "random start points; default: %u)" %
                          FbaParam.DEFAULT_NUMITER, metavar="N")
    optOptions.add_option("-l", "--use-full-matrix", action="store_true",
                          dest="useFullMatrix", help="use full matrix (disable "
                          "removal of dead ends, nonfunctional reactions, and "
                          "reactions with flux restricted to zero) - slow")
    parser.set_defaults(solver="default", useFullMatrix=False,
                        alpha=_DEFAULT_ALPHA, beta=_DEFAULT_BETA,
                        numIter=FbaParam.DEFAULT_NUMITER)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-w")
    parser.check_required("-o")

    if options.wMomaFvaFile:
        if options.alpha < 0.:
            print("Error: alpha must be non-negative.")
            exit()
        if options.beta <= 0.:
            print("Error: beta must be positive.")
            exit()

    if options.numIter < 1:
        print("Error: Number of NLP runs must be positive.")
        exit()

    # 2. Parse reaction file

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
        lb, ub = pparser.parse(options.paramFile)[-2:]
        linConstraints = pparser.lin_constraints
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

    # 4. Parse solution file for wildtype

    wtFlux = MetabolicFlux()
    try:
        wtFlux.readFromFile(options.wtSolution)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.wtSolution))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("An error occurred parsing file %s:" %
               os.path.basename(options.wtSolution))
        print(strerror)
        exit()

    if not wtFlux.hasSameReactions(model):
        print("Error: Solution and model must have the same reactions.")
        exit()

    # 5. (Optionally) Read FVA solution for reducing the solution space
    #    - skip if full matrix is to be used

    if options.redFvaFile and not options.useFullMatrix:
        try:
            fbaSolution, redMinmax = \
                FvAnalyzer.parseSolutionFile(options.redFvaFile)
        except IOError as strerror:
            print ("An error occurred while trying to read file %s:" %
                   os.path.basename(options.redFvaFile))
            print(strerror)
            exit()
        except SyntaxError as strerror:
            print ("Error in FVA solution file %s:" %
                   os.path.basename(options.redFvaFile))
            print(strerror)
            exit()

        if not fbaSolution.hasSameReactions(model):
            print ("Error in file %s: FVA solution and model must have the same"
                   " reactions." % os.path.basename(options.redFvaFile))
            exit()

        # Determine blocked reactions based on FVA solution
        blockedReactions = FvAnalyzer.getBlockedReactions(redMinmax)

    elif options.useFullMatrix:
        blockedReactions = []
    else:
        # Determine blocked reactions based on dead-end analysis
        blockedReactions = model.findDeadEnds(True)[1]

    # 6. (Optionally) Read FVA solution and compute weights for wMOMA

    if options.wMomaFvaFile:
        try:
            fbaSolution, wmomaMinmax = \
                FvAnalyzer.parseSolutionFile(options.wMomaFvaFile)
        except IOError as strerror:
            print ("An error occurred while trying to read file %s:" %
                   os.path.basename(options.wMomaFvaFile))
            print(strerror)
            exit()
        except SyntaxError as strerror:
            print ("Error in FVA solution file %s:" %
                   os.path.basename(options.wMomaFvaFile))
            print(strerror)
            exit()

        if not fbaSolution.hasSameReactions(model):
            print ("Error in file %s: FVA solution and model must have the same"
                   " reactions." % os.path.basename(options.wMomaFvaFile))
            exit()

        weights = MomaAnalyzer.getWeightsFromFluxVar(model, wmomaMinmax,
                                                     options.alpha, options.beta)

    # 7. Perform MOMA against the wildtype solution on the metabolic network

    moma = MomaAnalyzer(options.solver)
    if options.wMomaFvaFile:
        print("Performing wMOMA with alpha = %g, beta = %g." % (options.alpha,
                                                                options.beta))
        distance, solution, status, ndim = moma.runOnModel(model, wtFlux,
                                                           linConstraints, options.numIter, weights, blockedReactions)
    else:
        print("Performing MOMA...")
        distance, solution, status, ndim = moma.runOnModel(model, wtFlux,
                                                           linConstraints, options.numIter, blockedReactions=blockedReactions)

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print('\n'+'\n'.join(msgs))
    print ("Info: The reduced network has %u reactions and %u metabolites." %
           ndim[1::-1])
    print("Solver status:", printStatus(status))

    # 8. Write output to file

    if len(solution) != 0:
        print("\nOptimal value of objective function:", distance)
        diff = sum(solution.strictSqDiff(wtFlux).fluxDict.values())
        print("\nSquared distance to wildtype solution:", diff)
        if model.biomass_name:
            print("Biomass:", solution[model.biomass_name])
        print("Total absolute flux:", end=' ')
        print(abs(array(list(solution.fluxDict.values()))).sum())
        try:
            solution.writeToFile(options.outputFile)
        except IOError as strerror:
            print ("Unable to write to file %s:" %
                   os.path.basename(options.outputFile))
            print(strerror)
            exit()
    else:
        print("No output written.")
        exit()


if __name__ == "__main__":
    main()
