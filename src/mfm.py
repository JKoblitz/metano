#!/usr/bin/env python

""" metano Metabolite Flux Minimization

This script performs MFM on the FBA solution for a metabolic network (given as
a set of reaction equations), an objective function, and a number of inequality
constraints.

It first performs a standard FBA (or reads in an FBA solution) and then performs
variability analysis by the following algorithm:

1. Introduce a new constraint: Original objective function must have a value
                               of at least 95% of the optimal value (as
                               determined from the original FBA; threshold is
                               configurable).
2. Generate coefficient vectors expressing the (producing) metabolite fluxes as
   linear combinations of the reaction fluxes (flux variables).
3. Perform LP with the linear objective functions defined by these coefficient
   vectors: Successively minimize each objective function.


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

from defines import FbaParam, padNumber, COPYRIGHT_VERSION_STRING
from fba import OptionParser, FbAnalyzer
from reactionparser import ReactionParser
from paramparser import ParamParser
from metabolicmodel import MetabolicModel
from linearproblem import LinearProblem
from numpy import array, vstack, append, nan, isinf
import os

EPSILON = 1E-9
LARGE_NUM = 1E4

class MfmAnalyzer(object):
    """ Class for metabolite flux minimization
    """

    def __init__(self, solver="default"):
        """ initialize MFM class

        Keyword arguments:

        solver    -- LP solver to be used
        """
        self.solver = solver


    def run(self, objective, matrix, lb, ub, threshold=1.0, eqs=[], ineqs=[]):
        """ successively formulate and solve linear problems for each metabolite
            flux with inequality constraint 'objective value <= threshold'

        Arguments refer to split fluxes (i.e. non-negative flux variables).

        Keyword arguments:

        objective      -- coefficient vector for linear objective function
        matrix         -- stoichiometric matrix of the metabolic network
        lb             -- list of lower bounds (indexed like matrix columns)
        ub             -- list of upper bounds (indexed like matrix columns)
        threshold      -- objective function threshold
        eqs            -- list of (coefficient vector, right-hand side) pairs
                          for additional equality constraints
        ineqs          -- list of - '' - for additional inequality constraints

        Returns:  list of minimum values, indexed like metabolites
        """
        matrix = array(matrix)
        nMetabolites, nCols = matrix.shape
        if nMetabolites == 0:
            return []  # Nothing to do

        # Construct matrices of original equality and inequality constraints
        Aeq, beq, Aineq, bineq = FbAnalyzer.makeConstraintMatrices(matrix, eqs,
                                                                   ineqs)
        # Combine original inequality constraints and
        # 'objective function <= threshold'
        if Aineq is None:
            Aineq = [objective]
            bineq = [threshold]
        else:
            Aineq = vstack((Aineq, [objective]))
            bineq = append(bineq, threshold)

        result = []
        ubVec = []
        for i in range(nCols):
            # Impose artificial upper bounds so that all LPs are bounded
            if isinf(ub[i]):
                ubVec.append(LARGE_NUM)
            else:
                ubVec.append(ub[i])

        psSplit = LinearProblem(Aeq, beq, Aineq, bineq, array(lb), array(ubVec),
                                self.solver)

        # Compute coefficient vectors of producing metabolite fluxes
        # - pick only positive coefficients from stoichiometric matrix
        posMatrix = matrix.copy()
        for i in range(nMetabolites):
            for j in range(nCols):
                if posMatrix[i, j] < 0.:
                    posMatrix[i, j] = 0.

        # Successively minimize each flux
        for index in range(nMetabolites):
            try:
                psSplit.setObjective(map(float, posMatrix[index, :]))
            except Exception:
                # Skip all-zero rows
                result.append(0.)
                continue

            minval, s = psSplit.minimize()
            psSplit.resetObjective()

            if len(s) == 0:
                result.append(nan)
            else:
                result.append(minval)

        return result


    def runOnModel(self, model, fbaParams, threshp=.95, objFuncVal=None,
                   rmDeadEnds=True):
        """ perform metabolite flux minimization on the given model with
            objective value <= threshp * maximum

        Keyword arguments:

        model        -- the MetabolicModel
        fbaParams    -- FBA parameters
        threshp      -- threshold percentage (objective_value <= thresp*maximum)
        objFuncVal   -- FBA optimum for objective function (optional)
        rmDeadEnds   -- if True, remove all reactions with dead ends before
                        analysis (faster and gives an optimal solution, as well)

        Returns:
        minVec, dimReduced

        minVec       -- list of minimum values, indexed like metabolites
        dimReduced   -- pair (nRows, nColumns) with dimensions of reduced matrix
        """
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
                return [], array(modelRed.getStoichiometricMatrix()).shape
        else:
            modelRed = model

        matrix = array(modelRed.getStoichiometricMatrix())
        dimReduced = matrix.shape
        lb, ub = modelRed.getBounds()

        # Split fluxes into non-negative components
        matrixSplit, reactionsSplit, lbSplit, ubSplit = \
            FbAnalyzer.splitFluxes(matrix, modelRed.getReactionNames(), lb, ub)

        # Build (negative) objective function vector for split fluxes
        objective = ParamParser.convertObjFuncToLinVec(fbaParams.objStr,
            reactionsSplit, len(lbSplit), fbaParams.maxmin)
        maxmin_factor = -1. if fbaParams.maxmin else 1.

        # If the optimum of the objective function is not given, perform FBA
        if objFuncVal is None:
            fba = FbAnalyzer(self.solver)

            objFuncVal, sFlux = fba.run(reactionsSplit, matrixSplit, lbSplit,
                                        ubSplit, fbaParams)

            if len(sFlux) == 0:
                return [], dimReduced

        # Use negative threshold (objective value >= threshold is equivalent to
        # -objective value <= -threshold, and objective already has coefficient
        # -1 due to maximization)
        threshold = maxmin_factor * objFuncVal * threshp
        print "obj func. opt:", objFuncVal
        print "Threshold: ", threshold
        try:
            eqs, ineqs = ParamParser.linConstraintsToVectors(
                fbaParams.linConstraints, modelRed.reactionDict, len(lbSplit))
        except ValueError, e:
            # If any linear constraint is contradictory, report error
            print "Optimization not possible due to contradictory constraints:"
            print "  "+e
            exit()
        minVec = self.run(objective, matrixSplit, lbSplit, ubSplit, threshold,
            eqs, ineqs)


        # Add removed reactions with min = 0. to minVec and solution
        if len(modelRed) != len(model):
            minVecFull = []
            metabolitesRed = set(modelRed.getMetaboliteNames())

            for met in model.getMetaboliteNames():
                if met in metabolitesRed:
                    # Take values from minVec
                    index = modelRed.metaboliteDict[met]
                    minVecFull.append(minVec[index])
                else:
                    # Dead end => min = 0.
                    minVecFull.append(0.)

            return minVecFull, dimReduced
        else:
            return minVec, dimReduced


    @staticmethod
    def writeSolutionToFile(filename, metabolites, minVec):
        with open(filename, 'w') as f:
            MfmAnalyzer.writeSolutionToFileHandle(f, metabolites, minVec)

    @staticmethod
    def writeSolutionToFileHandle(f, metabolites, minVec):
        """ write MFM solution to the given file handle (must be open for
            writing)

        Keyword arguments:

        metabolites  -- dictionary { name : index }
        minVec       -- list of minimum flux values for each metabolite
        """
        if len(minVec) == 0:
            return  # Nothing to do
        if len(metabolites) != len(minVec):
            raise ValueError("Length mismatch: 'metabolites'(%u) 'minVec'(%u)"
                             % (len(metabolites), len(minVec)))

        # Construct output table and get widths of table columns (length of
        # longest occurring string)
        maxlenName, maxlenVal = 0, 0
        output = {}  #  dict {name : flux value as string }
        for met in metabolites:
            index = metabolites[met]
            # Pad non-negative numbers with space
            valStr = padNumber(repr(minVec[index]))
            output[met] = valStr
            lenName, lenVal = len(met), len(valStr)
            if lenName > maxlenName:
                maxlenName = lenName
            if lenVal > maxlenVal:
                maxlenVal = lenVal

        # Write table head
        f.write("NAME".ljust(maxlenName)+"   "+"MIN_FLUX".rjust(maxlenVal)+"\n")

        # Write result for the producing flux through every metabolite
        for met in sorted(metabolites):
            f.write(met.ljust(maxlenName)+" : "+output[met].rjust(maxlenVal)+
                    "\n")


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "Metabolite flux minimization\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform Metabolite Flux Minimization on the"
                      " network given by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use the given scenario FILE", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write MFM output to FILE", metavar="FILE")
    parser.add_option("-t", "--tolerance", dest="tolerance", type="float",
                      help="tolerance for objective function (0 < VALUE "
                      "< 1; default: .95)", metavar="VALUE")
    parser.add_option("-l", "--use-full-matrix", action="store_true",
                      dest="useFullMatrix", help="use full matrix (disable "
                      "removal of dead ends and nonfunctional reactions)")
    parser.set_defaults(tolerance=.95, useFullMatrix=False)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-o")

    if options.tolerance <= 0. or options.tolerance > 1.:
        print "Error: Tolerance must be in interval (0, 1]"
        exit()

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

    # 3. Parse scenario file

    model_messages = []
    pparser = ParamParser()
    try:
        # Parse file, get maxmin, name of objective function, and solver name
        maxmin, objStr, solver, numIter, lb , ub = pparser.parse(
                                                              options.paramFile)
        fbaParams = FbaParam(solver, maxmin, objStr, numIter)
        fbaParams.setLinConstraints(pparser.lin_constraints)
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
        print '\n'+'\n'.join(msgs)

    # 4. Launch metabolite flux minimization

    print "Launching metabolite flux minimization."
    print "tolerance: %g" % options.tolerance
    if solver == "":
        solver = "default"
    mfm = MfmAnalyzer(solver)
    minVec, ndim = mfm.runOnModel(model, fbaParams, options.tolerance,
                                  rmDeadEnds=not options.useFullMatrix)

    print ("Info: The reduced network has %u reactions and %u metabolites." %
           ndim[1::-1])
    if not minVec:
        print "Model is infeasible or unbounded. Nothing to do."
        exit()

    # 5. Write output to file

    try:
        mfm.writeSolutionToFile(options.outputFile, model.metaboliteDict,
                                minVec)
    except IOError, strerror:
        print ("Unable to write to file %s:" %
               os.path.basename(options.outputFile))
        print strerror
        exit()


if __name__ == "__main__":
    main()
