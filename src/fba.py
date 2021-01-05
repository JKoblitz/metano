#!/usr/bin/env python

""" metano Flux Balance Analysis

This script performs flux balance analysis (FBA) on a metabolic network
(given as a set of reaction equations), an objective function, and a number
of inequality constraints.

It parses a reaction file and a scenario file, formulates the
corresponding linear, quadratic, or other nonlinear optimization problem,
and calls an appropriate solver.


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

from builtins import range
from builtins import object
import optparse
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.linearproblem import LinearProblem, _GLPK, _OOGLPK
from metano.nonlinearproblem import NonLinearProblem
from metano.objectivefun import biomass_per_flux, neg_grad_biomass_per_flux
from metano.startpointiterator import StartPointIterator
from metano.defines import FbaParam, printHistogram, COPYRIGHT_VERSION_STRING
from numpy import array, insert, vstack, dot
import os


EPSILON = 1E-8


class OptionParser(optparse.OptionParser):
    """ Extension of OptionParser with function for checking required arguments
    """

    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            self.error("%s option not supplied" % option)


class FbAnalyzer(object):
    """ Class for flux-balance analysis
    """

    def __init__(self, solver="default", verbose=True):
        """ initialize FBA class

        Keyword arguments:

        solver    -- name of the solver to be used for the LP/QP/NLP
        """
        self.solver = solver
        self.verb = verbose

    @staticmethod
    def splitFluxes(matrix, reaction_names, lb, ub):
        """ split the fluxes into non-negative components

        This is done by duplicating and negating matrix columns for reversible
        reactions and just negating the matrix columns corresponding to flipped
        irreversible reactions.

        Keyword arguments:

        matrix         -- the stoichiometric matrix
        reaction_names -- list of reactions (indexed like matrix columns)
        lb             -- list of lower bounds (indexed like matrix columns)
        ub             -- list of upper bounds (indexed like matrix columns)

        Returns:
        matrixSplit, reactionsSplit, lbSplit, ubSplit

        matrixSplit    -- matrix with extra columns for negative flux components
        reactionsSplit -- dictionary { name : pair of indexes } - indexes can be
                          None; if both are not None, the first is the pos. and
                          the second is the neg. component (v = v[0] - v[1])
        lbSplit        -- list of lower bounds (indexed like columns of
                                                matrixSplit)
        ubSplit        -- list of upper bounds ( - '' - )
        """

        matrixSplit = matrix.copy()     # Perform deep copy
        reactionsSplit = {}
        lbSplit = []
        ubSplit = []
        newIndex = 0                    # Index in matrixSplit

        for i in range(len(reaction_names)):
            rea = reaction_names[i]
            if lb[i] >= 0.:
                # Reaction is irreversible (pos. flux only or restricted to 0)
                reactionsSplit[rea] = (newIndex, None)
                lbSplit.append(lb[i])
                ubSplit.append(ub[i])
            elif ub[i] <= 0.:
                # Reaction is irreversible & flipped => negate
                reactionsSplit[rea] = (None, newIndex)
                matrixSplit[:, newIndex] *= -1
                lbSplit.append(-ub[i])
                ubSplit.append(-lb[i])
            else:
                # Reaction is reversible => split (copy, then negate the copy)
                reactionsSplit[rea] = (newIndex, newIndex+1)
                matrixSplit = insert(matrixSplit, newIndex+1,
                                     -matrixSplit[:, newIndex], axis=1)
                lbSplit.append(0.)
                ubSplit.append(ub[i])
                lbSplit.append(0.)
                ubSplit.append(-lb[i])
                newIndex += 1

            newIndex += 1

        return matrixSplit, reactionsSplit, lbSplit, ubSplit

    @staticmethod
    def splitFluxVector(vec, reaction_names, reactionsSplit):
        """ split the given vector as described by reactionsSplit

        Keyword arguments:

        vec            -- the flux vector to be split
        reaction_names -- list of reactions corresponding to entries of vec
        reactionsSplit -- dictionary { name : tuple of indexes } - indexes can
                          be None; if both are not None, the first is the pos.
                          and the second is the neg. component (v = v[0] - v[1])
                          (as computed by FbAnalyzer.splitFluxes)

        Returns:  split vector (containing only non-negative values)
        """
        vecSplit = []
        for i in range(len(reaction_names)):
            rea = reaction_names[i]
            indexPos, indexNeg = reactionsSplit[rea]
            if len(reactionsSplit[rea]) != 2:
                raise TypeError("reactionsSplit tuples must exactly 2 entries")
            if vec[i] == 0.:
                if indexPos is not None:
                    vecSplit.append(0.)
                if indexNeg is not None:
                    vecSplit.append(0.)
            elif vec[i] > 0.:
                if indexPos is None:
                    print(rea, indexPos, indexNeg)
                    raise ValueError("Value %u (%g) out of range in split "
                                     "vector (must not be positive)" % (i, vec[i]))
                vecSplit.append(vec[i])
                if indexNeg is not None:
                    vecSplit.append(0.)
            else:
                if indexPos is not None:
                    vecSplit.append(0.)
                if indexNeg is None:
                    raise ValueError("Value %u (%g) out of range in split "
                                     "vector (must not be negative)" % (i, vec[i]))
                vecSplit.append(-vec[i])

        return vecSplit

    @staticmethod
    def rejoinFluxes(solutionSplit, reactionsSplit, reactions):
        """ rejoin the fluxes in the given solution by applying reactionsSplit

        This assumes that solutionSplit has been obtained by performing FBA on
        the split matrix corresponding to reactionsSplit. This function combines
        the split fluxes by subtracting the negative component from the positive
        one.

        Keyword arguments:

        solutionSplit  -- vector of non-negative fluxes
        reactionsSplit -- dictionary { name : tuple of indexes } - indexes can
                          be None; if both are not None, the first is the pos.
                          and the second is the neg. component (v = v[0] - v[1])
                          (as computed by FbAnalyzer.splitFluxes)
        reactions      -- dictionary of all reactions { name : original matrix
                                                                        column }

        Returns:  joined flux vector (indexed like original matrix columns)
        """
        if len(solutionSplit) == 0:
            return solutionSplit  # Nothing to do
        solution = [0.]*len(reactions)

        for rea in reactionsSplit:
            rea_index_orig = reactions[rea]
            indexPos, indexNeg = reactionsSplit[rea]

            if indexPos is not None:
                solution[rea_index_orig] += solutionSplit[indexPos]

            if indexNeg is not None:
                solution[rea_index_orig] -= solutionSplit[indexNeg]

        return solution

    @staticmethod
    def makeConstraintMatrices(matrix, eqs, ineqs):
        """ transform the given matrix and equality and inequality constraints
            to left-hand-side matrices and right-hand-side vectors of equality
            and inequality constraints

        Keyword arguments:

        matrix -- the stoichiometric matrix of the metabolic network
        eqs    -- list of (coefficient vector, right-hand side) pairs for the
                  equality constraints
        ineqs  -- list of - '' - for the inequality constraints

        """
        if eqs:
            Aeq = vstack((matrix, array([row[0] for row in eqs])))
            beq = array([0.]*len(matrix) + [row[1] for row in eqs])
        else:
            Aeq = matrix
            beq = array([0.]*len(matrix))
        if ineqs:
            Aineq = array([row[0] for row in ineqs])
            bineq = array([row[1] for row in ineqs])
        else:
            Aineq, bineq = None, None

        return Aeq, beq, Aineq, bineq

    def run(self, reactions, matrix, lb, ub, fbaParams):
        """ analyze objective function, and solve the LP/QP/NLP

        Keyword arguments:

        reactions  -- dictionary of all reactions { name : matrix column } or
                      { name : tuple of indices } for split fluxes
        matrix    -- the stoichiometric matrix of the metabolic network
        lb        -- list of lower bounds (indexed like matrix columns)
        ub        -- list of upper bounds (indexed like matrix columns)
        fbaParams -- optimization parameters (incl. objective function)

        Returns:
        obj_value, solution

        obj_value -- optimal value of objective function
        solution  -- a solution where the objective function assumes obj_value
        """
        maxmin, objStr, numIter = (fbaParams.maxmin, fbaParams.objStr,
                                   fbaParams.numIter)
        # Multiply by -1 for maximization
        maxmin_factor = -1. if maxmin == True else 1.
        nCols = len(lb)

        # Get additional linear equality and inequality constraints
        try:
            Aeq, beq, Aineq, bineq = self.makeConstraintMatrices(matrix,
                                                                 *ParamParser.linConstraintsToVectors(fbaParams.linConstraints,
                                                                                                      reactions, nCols))
        except ValueError:
            # If any linear constraint is contradictory, return empty solution
            return 0., []

        if objStr.lower().startswith("per_flux"):
            # special nonlinear objective function: single flux per flux unit
            # e.g. "per_flux(Biomass)"
            try:
                str_perflux, str_flux = objStr.split('(', 2)
            except ValueError:
                print ("Error in scenario file in objective function "
                       "definition.\nPER_FLUX syntax: "
                       "PER_FLUX(<reaction_name>)")
                exit()

            if str_perflux.strip().lower() != "per_flux":
                print ("Error: Objective function '%s' is not a reaction name "
                       "or a PER_FLUX definition.\n"
                       "Currently, objective function must be a linear "
                       "combination of reaction fluxes or PER_FLUX(<reaction>)."
                       % objStr)
                exit()

            rea_name = str_flux.split(')', 1)[0].strip()
            try:
                biomass_index = reactions[rea_name][0]
            except TypeError:
                biomass_index = reactions[rea_name]

            def obj_func(x): return -biomass_per_flux(x, biomass_index)

            if numIter < 0:
                numIter = FbaParam.DEFAULT_NUMITER

            try:
                ps = NonLinearProblem(Aeq, beq, Aineq, bineq, lb, ub,
                                      self.solver)
                ps.setObjective(obj_func)
                ps.setObjGrad(lambda x:
                              neg_grad_biomass_per_flux(x, biomass_index))
                if fbaParams.nlc != []:
                    ps.setNonlinearConstraints(fbaParams.nlc)
                    ps.setNlcGrad(fbaParams.nlc_grad)
                spIter = StartPointIterator(nCols, numIter)
                spIter.setRange(-1., 1.)
                ps.setStartPointIterator(spIter)

            except ValueError as strerror:
                print(strerror)
                exit()

        else:
            # Build linear objective function (given as coefficient vector)
            try:
                objective = ParamParser.convertObjFuncToLinVec(objStr,
                                                               reactions, nCols, maxmin)
            except Exception as strerror:
                print ("Error while trying to build coefficient vector of "
                       "linear objective function:")
                print(strerror)
                exit()

            if fbaParams.nlc == []:
                # Linear function and linear constraints only

                # If flux variables are not split, use GLPK through OpenOpt
                solver = (_OOGLPK if self.solver.lower() == _GLPK and
                          nCols == len(reactions) else self.solver)
                try:
                    ps = LinearProblem(Aeq, beq, Aineq, bineq, lb, ub, solver)
                except ValueError as strerror:
                    print(strerror)
                    exit()
                #print "lb:", len(lb)
                ps.setObjective(objective)
            else:
                # Linear function but nonlinear constraints
                try:
                    ps = NonLinearProblem(Aeq, beq, Aineq, bineq, lb, ub,
                                          self.solver)
                except ValueError as strerror:
                    print(strerror)
                    exit()

                ps.setObjective(lambda x: dot(x, objective))
                ps.setObjGrad(lambda _: array(objective))
                ps.setNonlinearConstraints(fbaParams.nlc)
                ps.setNlcGrad(fbaParams.nlc_grad)
                spIter = StartPointIterator(nCols, numIter)
                spIter.setRange(-1., 1.)
                ps.setStartPointIterator(spIter)

        # Launch optimization problem solver
        try:
            obj_value, solution = ps.minimize()
            print(obj_value)
        except ValueError as strerror:
            if self.solver == "default":
                print("Default solver reported an error:")
            else:
                print("Solver %s reported an error:" % self.solver)
            print(strerror)
            print("Consider trying a different solver.")
            exit()
        except KeyError as strerror:
            print(strerror)
            exit()

        obj_value *= maxmin_factor
        if hasattr(ps, "getHistogram"):
            printHistogram(*ps.getHistogram(32)+(maxmin_factor < 0.,))

        return obj_value, solution

    def runWithTotalFluxMinimization(self, reactions, matrix, lb, ub, fbaParams,
                                     start_value_total=None, tolerance_total=EPSILON, tolerance_obj=EPSILON):
        """ perform FBA with LP iteratively to find a flux distribution with
            optimal value of objective function value and minimum total flux

        This only works with non-negative flux variables, i.e. split fluxes!
        The objective function and all constraints must be linear.

        This function limits the total flux to a parameter and then iteratively
        fits this parameter until either the change in total flux is less than
        tolerance_total or the change in the objective function is less than
        tolerance_obj

        Keyword arguments:

        reactions  -- dictionary { name : tuple of indices } for split fluxes
        matrix    -- the stoichiometric matrix of the metabolic network
        lb        -- list of lower bounds (indexed like matrix columns)
        ub        -- list of upper bounds (indexed like matrix columns)
        fbaParams -- optimization parameters (incl. objective function)
        start_value_total -- initial flux limit (if None: set to len(reactions))
        tolerance_total   -- stopping criterion: change in flux limit below this
        tolerance_obj     -- stopping criterion: change in objective function
                                                 value below this

        Returns:
        limit, obj_value, solution, nSteps

        limit     -- limit of total flux
        obj_value -- optimal value of objective function
        solution  -- a solution where the objective function assumes obj_value,
                     obtained with limited total flux
        nSteps    -- number of steps until optimal limit was found
        """
        maxmin, objStr = (fbaParams.maxmin, fbaParams.objStr)
        # Multiply by -1 for maximization
        maxmin_factor = -1. if maxmin == True else 1.
        nCols = len(lb)
        for value in lb:
            if value < 0.:
                print ("Error: Minimization of total flux requires non-negative"
                       " flux variables.")
                exit()

        # Get additional linear equality and inequality constraints
        try:
            Aeq, beq, Aineq, bineq = self.makeConstraintMatrices(matrix,
                                                                 *ParamParser.linConstraintsToVectors(fbaParams.linConstraints,
                                                                                                      reactions, nCols))
        except ValueError:
            # If any linear constraint is contradictory, return empty solution
            return 0., []

        # Prepare additional constraint: total flux < threshold
        if Aineq is None:
            Aineq, bineq = array([[1.]*nCols]), []
        else:
            # Aineq has extra line for total flux
            Aineq = vstack((Aineq, [[1.]*nCols]))
            # Convert bineq to list (allows easy addition of extra value)
            bineq = list(bineq)

        # Build linear objective function (given as coefficient vector)
        try:
            objective = ParamParser.convertObjFuncToLinVec(objStr, reactions,
                                                           nCols, maxmin)
        except Exception as strerror:
            print ("Error while trying to build coefficient vector of "
                   "linear objective function:")
            print(strerror)
            exit()

        if start_value_total is None or start_value_total < 0.:
            start_value_total = float(len(reactions))

        # Prepare left and right border for binary search
        l = start_value_total / 2.
        r = start_value_total * 1.5
        m = start_value_total

        # Use GLPK through OpenOpt
        solver = (_OOGLPK if (self.solver.lower() == _GLPK or
                              self.solver.lower() == "default") else self.solver)
        try:
            ps = LinearProblem(Aeq, beq, Aineq, bineq+[m], lb, ub, solver)
        except ValueError as strerror:
            print(strerror)
            exit()
        ps.setObjective(objective)

        try:
            m_value = ps.minimize()[0]
        except ValueError as strerror:
            if self.solver == "default":
                print("Default solver reported an error:")
            else:
                print("Solver %s reported an error:" % self.solver)
            print(strerror)
            print("Consider trying a different solver.")
            exit()
        except KeyError as strerror:
            print(strerror)
            exit()

        ps.bineq[-1] = l
        l_value, l_flux = ps.minimize()
        #print "l_flux", l_flux
        ps.bineq[-1] = r
        r_value = ps.minimize()[0]

        nSteps = 3     # Counter for calls of ps.minimize()
        print ("Step %u: l=%.10g (%.10g), r=%.10g (%.10g), m=%.10g (%.10g)" %
               (nSteps, l, l_value, r, r_value, m, m_value))

        # 1. Find right border so that f(m) = f(r)
        while abs(m_value-r_value) >= tolerance_obj:
            r *= 10.
            ps.bineq[-1] = r
            r_value = ps.minimize()[0]
            m = (l+r)/2.
            ps.bineq[-1] = m
            m_value = ps.minimize()[0]
            nSteps += 2
            print ("Step %u: l=%.10g (%.10g), r=%.10g (%.10g), m=%.10g (%.10g)"
                   % (nSteps, l, l_value, r, r_value, m, m_value))
        print("Right border found.")
        print("abs(m_value-r_value):", abs(m_value-r_value))

        # 2. Find left border so that f(l) > f(m)
        l_changed = False
        while l_flux == [] or (abs(l_value-m_value) < tolerance_obj and
                               abs(l_value) >= EPSILON):
            if not l_changed:
                l_changed = True
            if l_flux == []:
                l = (l+m)/2.
            else:
                l /= 2.
            ps.bineq[-1] = l
            l_value, l_flux = ps.minimize()
            nSteps += 1

        if l_changed:
            m = (l+r)/2.
            ps.bineq[-1] = m
            m_value = ps.minimize()[0]
            nSteps += 1
            print ("Step %u: l=%.10g (%.10g), r=%.10g (%.10g), m=%.10g (%.10g)"
                   % (nSteps, l, l_value, r, r_value, m, m_value))
        print("Left border found.")
        print("abs(m_value-l_value):", abs(m_value-l_value))

        while (abs(r-l) >= tolerance_total and
               abs(l_value-r_value) >= tolerance_obj):

            if abs(m_value-r_value) < tolerance_obj:
                r, r_value = m, m_value

            else:
                l, l_value = m, m_value

            m = (l+r)/2.
            ps.bineq[-1] = m
            m_value, solution = ps.minimize()
            nSteps += 1

            print ("Step %u: l=%.10g (%.10g), r=%.10g (%.10g), m=%.10g (%.10g)"
                   % (nSteps, l, l_value, r, r_value, m, m_value))

        if r_value < m_value:
            m = r
            ps.bineq[-1] = m
            m_value, solution = ps.minimize()
            nSteps += 1
            print ("Step %u: l=%.10g (%.10g), r=%.10g (%.10g), m=%.10g (%.10g)"
                   % (nSteps, l, l_value, r, r_value, m, m_value))

        return m, m_value * maxmin_factor, solution, nSteps

    def runOnModel(self, model, fbaParams, splitFluxes=True,
                   minimizeTotalFlux=False, rmDeadEnds=True):
        """ run flux-balance analysis on the given MetabolicModel

        Keyword arguments:

        model       -- the MetabolicModel with reactions and bounds
        fbaParams   -- optimization parameters (incl. objective function)
        splitFluxes -- if True, run split fluxes (resulting in non-negative flux
                       variables) before analysis
        minimizeTotalFlux -- if True, find a solution with minimal total absolute flux
        rmDeadEnds   -- if True, remove all reactions with dead ends before
                        analysis (faster and gives an optimal solution, as well)

        Returns:
        obj_value, solution, dimReduced

        obj_value   -- optimal value of objective function
        solution    -- a solution where the objective function assumes obj_value
        dimReduced  -- pair (nRows, nColumns) with dimensions of reduced matrix
        """
        maxmin = fbaParams.maxmin
        objStr = fbaParams.objStr

        if isinstance(maxmin, str):
            maxmin = {"max": True, "min": False}[maxmin.lower()]

        if rmDeadEnds:
            deadReactions = model.findDeadEnds(True)[1]
            modelRed = model.getSubModelByExcludeList(deadReactions)
            cbz = model.canBeZero(deadReactions)
            nonZeroDeadEnds = [deadReactions[i] for i in
                               range(len(deadReactions)) if not cbz[i]]
            if nonZeroDeadEnds:
                print ("The following blocked reactions are constrained to a "
                       "non-zero flux:\n  " + "\n  ".join(nonZeroDeadEnds) +
                       "\nThe problem is infeasible.")
                return (0., MetabolicFlux(),
                        array(modelRed.getStoichiometricMatrix()).shape)
        else:
            modelRed = model

        matrix = array(modelRed.getStoichiometricMatrix())
        dimReduced = matrix.shape
        reaction_names = modelRed.getReactionNames()
        reactions = modelRed.reactionDict
        lb, ub = modelRed.getBounds()

        if splitFluxes:
            matrixSplit, reactionsSplit, lbSplit, ubSplit = \
                self.splitFluxes(matrix, reaction_names, lb, ub)
            try:
                ParamParser.convertObjFuncToLinVec(objStr, reactionsSplit,
                                                   len(lbSplit), maxmin)
            except KeyError as s:
                if deadReactions:
                    if self.verb:
                        print ("Error while trying to construct the objective "
                               "function vector:\n%s\nThis may be due to removal "
                               "of nonfunctional reactions." % s)
                    return 0., MetabolicFlux(), dimReduced

            obj_value, solutionSplit = self.run(reactionsSplit, matrixSplit,
                                                lbSplit, ubSplit, fbaParams)
            solution = []
            try:
                if solutionSplit == []:
                    solution = []
                else:
                    if minimizeTotalFlux:
                        # Optimize solution by limiting total flux
                        limit, obj_value, solutionSplit, nSteps = \
                            self.runWithTotalFluxMinimization(reactionsSplit,
                                                              matrixSplit, lbSplit, ubSplit, fbaParams,
                                                              sum(solutionSplit), 1e-8, 1e-12)
                        print ("Found optimal flux limit of %.10g in %u steps."
                               % (limit, nSteps))

                    solution = array(FbAnalyzer.rejoinFluxes(solutionSplit,
                                                             reactionsSplit, reactions))
            except NameError as strerror:
                print(strerror)

        else:
            obj_value, solution = self.run(
                reactions, matrix, lb, ub, fbaParams)
        flux = MetabolicFlux(modelRed, solution)

        # Add removed reactions with flux 0. and original bounds to solution
        if len(modelRed) != len(model) and len(solution) != 0:
            reactionsRed = set(reaction_names)
            for rea in model:
                if rea.name not in reactionsRed:
                    flux.fluxDict[rea.name] = 0.
                    flux.boundsDict[rea.name] = (rea.lb, rea.ub)

        return obj_value, flux, dimReduced


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "Flux balance analysis\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform Flux Balance Analysis on the network given "
                           "by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use the given scenario FILE for Flux Balance "
                           "Analysis", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write output of Flux Balance Analysis to FILE",
                      metavar="FILE")
    parser.add_option("-m", "--minimize-total", action="store_true",
                      dest="minimizeFlux", help="secondarily minimize total "
                      "flux")
    parser.add_option("-l", "--use-full-matrix", action="store_true",
                      dest="useFullMatrix", help="use full matrix (disable "
                      "removal of dead ends and nonfunctional reactions)")
    parser.set_defaults(minimizeFlux=False, useFullMatrix=False)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-o")

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
        # Parse file, get maxmin, name of objective function, and solver name
        maxmin, objStr, solver, numIter, lb, ub = pparser.parse(
            options.paramFile)
        fbaParams = FbaParam(solver, maxmin, objStr, numIter)
        fbaParams.setLinConstraints(pparser.lin_constraints)
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

    # 4. Perform flux balance analysis on the metabolic network

    if solver == "":
        solver = "default"
    fba = FbAnalyzer(solver)
    obj_value, solution, ndim = fba.runOnModel(model, fbaParams, True,
                                               options.minimizeFlux, not options.useFullMatrix)

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print('\n'+'\n'.join(msgs))
    print ("Info: The reduced network has %u reactions and %u metabolites." %
           ndim[1::-1])

    # 5. Write output to file
    if len(solution) != 0:
        print("\nValue of objective function at solution:", obj_value)
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
