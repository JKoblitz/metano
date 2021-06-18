#!/usr/bin/env python

""" metano Flux Variability Analyzer

This script performs flux variability analysis on the FBA solution for a
metabolic network (given as a set of reaction equations), an objective
function, and a number of inequality constraints.

It first performs a standard FBA (or reads in an FBA solution) and then performs
variability analysis by the following algorithm:

- Introduce a new constraint: Original objective function must have a value
                              of at least 95% of the optimal value (as
                              determined from the original FBA; threshold is
                              configurable)
- Perform LP for each flux variable: Successively maximize and minimize each of
  the flux variables


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
import os
import csv
from metano.defines import FbaParam, padNumber, COPYRIGHT_VERSION_STRING
from metano.fba import OptionParser, FbAnalyzer
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.linearproblem import LinearProblem
from numpy import array, vstack, append, dot, nan, isinf


EPSILON = 1E-9
LARGE_NUM = 1E5


class FvAnalyzer(object):
    """ Class for flux variability analysis
    """

    def __init__(self, solver="default"):
        """ initialize FVA class

        Keyword arguments:

        solver    -- LP solver to be used
        """
        self.solver = solver

    def run(self, objective, matrix, lb, ub, threshold=.95, eqs=[], ineqs=[]):
        """ successively formulate and solve linear problems for each flux with
            inequality constraint 'objective value <= threshold'

        Keyword arguments:

        objective -- coefficient vector for linear objective function
        matrix    -- stoichiometric matrix
        lb        -- list of lower bounds, indexed like reactions
        ub        -- list of upper bounds, indexed like reactions
        threshold -- objective function threshold
        eqs       -- list of (coefficient vector, right-hand side) pairs for
                     additional equality constraints
        ineqs     -- list of - '' - for additional inequality constraints

        Returns:  list of pairs (minimum, maximum), indexed like reactions
        """
        # Construct matrices of original equality and inequality constraints
        Aeq, beq, Aineq, bineq = FbAnalyzer.makeConstraintMatrices(matrix, eqs,
                                                                   ineqs)
        #lb = [-1000]*939
        #ub = [1000]*939
        # Combine original inequality constraints and
        # 'objective function <= threshold'
        if Aineq is None:
            Aineq = [objective]
            bineq = [threshold]
        else:
            Aineq = vstack((Aineq, [objective]))
            bineq = append(bineq, threshold)

        nReactions = len(lb)
        lbVec, ubVec = [], []
        for i in range(nReactions):
            if isinf(lb[i]):
                lbVec.append(-LARGE_NUM)
            else:
                lbVec.append(lb[i])
            if isinf(ub[i]):
                ubVec.append(LARGE_NUM)
            else:
                ubVec.append(ub[i])

        ps = LinearProblem(Aeq, beq, Aineq, bineq, lbVec, ubVec, self.solver)
        minmax = [None]*nReactions

        # First maximize each flux
        for index in range(nReactions):

            ps.setObjective([0.]*index+[-1.]+[0.]*(nReactions-index-1))
            s = ps.minimize()[1]
            ps.resetObjective()

            # catch random infeasible solution that resolves if the solver is recreated
            if len(s) == 0:
                ps = LinearProblem(Aeq, beq, Aineq, bineq,
                                   lbVec, ubVec, self.solver)
                ps.setObjective([0.]*index+[-1.]+[0.]*(nReactions-index-1))
                s = ps.minimize()[1]
                ps.resetObjective()

            if len(s) != 0:
                maxval = s[index]
            else:
                maxval = nan
                s = None
            minmax[index] = maxval

        # Now minimize each flux
        for index in range(nReactions):

            ps.setObjective([0.]*index+[1.]+[0.]*(nReactions-index-1))
            s = ps.minimize()[1]
            ps.resetObjective()

            # catch random infeasible solution that resolves if the solver is recreated
            if len(s) == 0:
                ps = LinearProblem(Aeq, beq, Aineq, bineq,
                                   lbVec, ubVec, self.solver)
                ps.setObjective([0.]*index+[-1.]+[0.]*(nReactions-index-1))
                s = ps.minimize()[1]
                ps.resetObjective()

            if len(s) != 0:
                minval = s[index]
            else:
                minval = nan
                s = None
            minmax[index] = minval, minmax[index]

        return minmax

    def runOnModel(self, model, fbaParams, threshp=.95, solution=None,
                   splitFluxes=True, rmDeadEnds=True):
        """ perform flux variability analysis on the given model with
            objective value <= threshp * maximum

        Keyword arguments:

        model        -- the MetabolicModel
        fbaParams    -- FBA parameters
        threshp      -- threshold percentage (objective_value <= thresp*maximum)
        solution     -- optional: either FBA solution file or solution as
                        MetabolicFlux object
        splitFluxes  -- if True, run split fluxes (resulting in non-negative
                        flux variables) before FBA
        rmDeadEnds   -- if True, remove all reactions with dead ends before
                        analysis (faster and gives an optimal solution, as well)

        Returns:
        minmax, sFlux, lb, ub, dimReduced

        minmax       -- list of pairs (minimum, maximum), indexed like reactions
        sFlux        -- FBA solution (vector indexed like reactions)
        lb, ub       -- lower/upper bounds vectors, indexed like reactions
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
                lbFull, ubFull = list(map(array, model.getBounds()))
                return ([], MetabolicFlux(), lbFull, ubFull,
                        array(modelRed.getStoichiometricMatrix()).shape)
        else:
            modelRed = model

        matrix = array(modelRed.getStoichiometricMatrix())
        dimReduced = matrix.shape
        lb, ub = list(map(array, modelRed.getBounds()))

        # Build (negative) objective function vector
        objective = ParamParser.convertObjFuncToLinVec(fbaParams.objStr,
                                                       modelRed.reactionDict, -1, fbaParams.maxmin)
        maxmin_factor = -1. if fbaParams.maxmin else 1.

        # Case 1. If an FBA solution is given, use that

        if isinstance(solution, MetabolicFlux):
            # Sort solution like matrix columns
            sFlux = solution.getVecOrderedByModel(modelRed)
            # Evaluate the objective function at solution
            obj_value = maxmin_factor*dot(objective, sFlux)

        # Case 2. If solution file is given, read solution

        elif solution is not None:
            sFlux = MetabolicFlux()
            try:
                sFlux.readFromFile(solution)
            except IOError as strerror:
                print ("An error occurred while trying to read file %s:" %
                       os.path.basename(solution))
                print(strerror)
                exit()
            except SyntaxError as strerror:
                print ("An error occurred parsing file %s:" %
                       os.path.basename(solution))
                print(strerror)
                exit()

            if not sFlux.hasSameReactions(model):
                print("Error: Solution and model must have the same reactions.")
                exit()

            # Sort solution like matrix columns
            sFlux = sFlux.getVecOrderedByModel(modelRed)

            # Evaluate the objective function at solution
            obj_value = maxmin_factor*dot(objective, sFlux)

        # Case 3. Perform flux balance analysis on the metabolic network
        #         (only if no solution file is given)

        else:
            fba = FbAnalyzer(self.solver)

            if splitFluxes:
                # Split flux variables for FBA
                matrixSplit, reactionsSplit, lbSplit, ubSplit = \
                    FbAnalyzer.splitFluxes(matrix, modelRed.getReactionNames(),
                                           lb, ub)
                obj_value, sFlux = fba.run(reactionsSplit, matrixSplit, lbSplit,
                                           ubSplit, fbaParams)
                # Get joined solution (for output)
                sFlux = array(FbAnalyzer.rejoinFluxes(sFlux, reactionsSplit,
                                                      modelRed.reactionDict))
            else:
                obj_value, sFlux = fba.run(modelRed.getReactionNames(), matrix,
                                           lb, ub, fbaParams)

            if len(sFlux) == 0:
                lbFull, ubFull = list(map(array, model.getBounds()))
                return [], sFlux, lbFull, ubFull, dimReduced

        # Use negative threshold (objective value >= threshold is equivalent to
        # -objective value <= -threshold, and objective already has coefficient
        # -1 due to maximization)
        threshold = maxmin_factor * obj_value * threshp
        print("obj func. opt:", obj_value)
        print("Threshold: ", threshold)
        try:
            eqs, ineqs = ParamParser.linConstraintsToVectors(
                fbaParams.linConstraints, modelRed.reactionDict)
        except ValueError as e:
            # If any linear constraint is contradictory, report error
            print("Optimization not possible due to contradictory constraints:")
            print("  "+e)
            exit()

        minmax = self.run(objective, matrix, lb, ub, threshold, eqs, ineqs)

        # Add removed reactions with min = max = 0. to minmax and solution
        if len(modelRed) != len(model):
            minmaxFull, sFluxFull = [], []
            reactionsRed = set(modelRed.getReactionNames())

            for rea in model:
                if rea.name in reactionsRed:
                    # Take values from minmax and sFlux
                    index = modelRed.reactionDict[rea.name]
                    minmaxFull.append(minmax[index])
                    sFluxFull.append(sFlux[index])
                else:
                    # Dead reaction => min = max = flux = 0.
                    minmaxFull.append((0., 0.))
                    sFluxFull.append(0.)

            lbFull, ubFull = list(map(array, model.getBounds()))
            return minmaxFull, array(sFluxFull), lbFull, ubFull, dimReduced
        else:
            return minmax, sFlux, lb, ub, dimReduced

    @staticmethod
    def writeSolutionToFile(filename, reactions, solution, lb, ub, minmax):
        with open(filename, 'w') as f:
            if filename.split(".")[-1] == "csv":
                FvAnalyzer.writeSolutionToCSVHandle(f, reactions, solution, lb, ub,
                                                    minmax)
            else:
                FvAnalyzer.writeSolutionToFileHandle(f, reactions, solution, lb, ub,
                                                     minmax)

    @staticmethod
    def writeSolutionToCSVHandle(f, reactions, solution, lb, ub, minmax):
        """ write to the CSV file given by file handle f (must be open for writing)
        """
        if len(solution) == 0:
            return  # Nothing to do

        csvwriter = csv.writer(
            f, delimiter=";", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        # Write Header to CSV file
        csvwriter.writerow(["NAME", "MIN_FLUX", "FBA_FLUX",
                            "MAX_FLUX", "DIFF", "LB", "UB"])
        # Construct output table and get widths of table columns (length of
        # longest occurring string)
        for rea in sorted(reactions):
            index = reactions[rea]
            mini, maxi = minmax[index]
            diff = maxi-mini
            flux = solution[index]
            lbVec = lb[index]
            ubVec = ub[index]

            csvwriter.writerow([rea, mini, flux, maxi, diff, lbVec, ubVec])

    @staticmethod
    def writeSolutionToFileHandle(f, reactions, solution, lb, ub, minmax):
        if len(solution) == 0:
            return  # Nothing to do

        # Construct output table and get widths of table columns (length of
        # longest occurring string)
        maxlenName, maxlenMin, maxlenFlux, maxlenMax, maxlenDiff = 0, 0, 0, 0, 0
        maxlenLb, maxlenUb = 0, 0
        #  dict {name : (minFlux, fbaFlux, maxFlux, diff, lb, ub)} - as strings
        outputTable = {}
        for rea in reactions:
            index = reactions[rea]
            mini, maxi = minmax[index]
            # Pad non-negative numbers with space
            outputTable[rea] = tuple(map(padNumber, list(map(repr, (solution[index],
                                                                    lb[index], ub[index], mini, maxi, maxi-mini)))))
            lenName, lenFlux, lenLb, lenUb, lenMin, lenMax, lenDiff = \
                list(map(len, (rea,)+outputTable[rea]))

            if lenName > maxlenName:
                maxlenName = lenName
            if lenFlux > maxlenFlux:
                maxlenFlux = lenFlux
            if lenLb > maxlenLb:
                maxlenLb = lenLb
            if lenUb > maxlenUb:
                maxlenUb = lenUb
            if lenMin > maxlenMin:
                maxlenMin = lenMin
            if lenMax > maxlenMax:
                maxlenMax = lenMax
            if lenDiff > maxlenDiff:
                maxlenDiff = lenDiff

        # Write table head
        f.write("NAME".ljust(maxlenName)+"   "+"MIN_FLUX".rjust(maxlenMin)+" " +
                "FBA_FLUX".rjust(maxlenFlux)+" " +
                "MAX_FLUX".rjust(maxlenMax)+" "
                + "DIFF".rjust(maxlenDiff)+" "+"LB".rjust(maxlenLb)+" " +
                "UB".rjust(maxlenUb)+"\n")

        # Write result for the flux through every reaction
        for rea in sorted(reactions):
            flux, lb, ub, mini, maxi, diff = outputTable[rea]
            f.write(rea.ljust(maxlenName) + " : " + mini.rjust(maxlenMin) + " "
                    + flux.rjust(maxlenFlux) + " " +
                    maxi.rjust(maxlenMax) + " "
                    + diff.rjust(maxlenDiff) + " " + lb.rjust(maxlenLb) + " " +
                    ub.rjust(maxlenUb) + "\n")

    @staticmethod
    def parseSolutionFile(filename):
        """ read the given FVA solution file identified by filename
            -- wrapper for parseSolutionFileByHandle()
        """
        with open(filename) as f:
            if filename.split(".")[-1] == "csv":
                return FvAnalyzer.parseCSVFileByHandle(f)
            else:
                return FvAnalyzer.parseSolutionFileByHandle(f)

    @staticmethod
    def parseCSVFileByHandle(f):
        """ read the FVA solution file given as a file object (must be open)

        Returns:
        (fbaSolution, minmax) tuple with

        fbaSolution -- MetabolicFlux object with flux distribution, and lb, ub
        minmax      -- dict { reaction : flux minimum, flux maximum }
        """
        csvreader = csv.reader(f, delimiter=";", quotechar='"')
        reactions = set()
        reactionVec, fluxVec = [], []
        boundsDict = {}
        minmax = {}

        line_no = 0
        for line in csvreader:
            line_no += 1
            if line[0] == "NAME":
                continue

            rea = line[0]
            if rea in reactions:
                raise SyntaxError("Syntax error in line %u: Duplicate "
                                  "reaction." % line_no)
            reactions.add(rea)
            reactionVec.append(rea)

            try:
                values = list(map(float, line[1:]))
            except ValueError:
                raise SyntaxError("Syntax error in line %u: Invalid "
                                  "floating point value." % line_no)
            try:
                fluxVec.append(values[1])
                minmax[rea] = values[0], values[2]
                boundsDict[rea] = values[4], values[5]
            except IndexError:
                raise SyntaxError("Syntax error in line %u:\nLine must "
                                  "contain exactly six values (min_flux, "
                                  "fba_flux, max_flux, diff, lb, ub)." %
                                  line_no)
        return MetabolicFlux(reactionVec, fluxVec, boundsDict), minmax

    @staticmethod
    def parseSolutionFileByHandle(f):
        """ read the FVA solution file given as a file object (must be open)

        Returns:
        (fbaSolution, minmax) tuple with

        fbaSolution -- MetabolicFlux object with flux distribution, and lb, ub
        minmax      -- dict { reaction : flux minimum, flux maximum }
        """
        reactions = set()
        reactionVec, fluxVec = [], []
        boundsDict = {}
        minmax = {}

        line_no = 0
        for line in f:
            line_no += 1
            if (line.lstrip().upper().startswith("NAME") or line == "" or
                    line.isspace()):
                continue

            try:
                rea, values_str = list(map(str.rstrip, line.split(":")))
            except ValueError:
                raise SyntaxError("Syntax error in line %u:\nLine must "
                                  "contain exactly one colon (':')." %
                                  line_no)
            if rea in reactions:
                raise SyntaxError("Syntax error in line %u: Duplicate "
                                  "reaction." % line_no)
            reactions.add(rea)
            reactionVec.append(rea)

            try:
                values = list(map(float, values_str.split(None, 6)[:6]))
            except ValueError:
                raise SyntaxError("Syntax error in line %u: Invalid "
                                  "floating point value." % line_no)
            try:
                fluxVec.append(values[1])
                minmax[rea] = values[0], values[2]
                boundsDict[rea] = values[4], values[5]
            except IndexError:
                raise SyntaxError("Syntax error in line %u:\nLine must "
                                  "contain exactly six values (min_flux, "
                                  "fba_flux, max_flux, diff, lb, ub)." %
                                  line_no)
        return MetabolicFlux(reactionVec, fluxVec, boundsDict), minmax

    @staticmethod
    def getBlockedReactions(fvaMinmax):
        """ get reactions whose flux is restricted to zero based on FVA solution

        Keyword arguments:
        fvaMinmax  -- dict { reaction : flux minimum, flux maximum }

        Returns list of reactions (in arbitrary order)
        """
        blockedReactions = []
        for name in fvaMinmax:
            mini, maxi = fvaMinmax[name]
            if mini == maxi == 0.:
                blockedReactions.append(name)
        return blockedReactions


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "Flux variability analysis\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform Flux Variability Analysis on the network "
                           "given by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use the given scenario FILE for Flux Variability "
                           "Analysis", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write output of Flux Variability Analysis to FILE",
                      metavar="FILE")
    parser.add_option("-t", "--tolerance", dest="tolerance", type="float",
                      help="tolerance for objective function (0 < VALUE "
                      "< 1; default: .95)", metavar="VALUE")
    parser.add_option("-s", "--solution-file", dest="solutionFile", help="read "
                      "FBA solution from FILE (optional)", metavar="FILE")
    parser.add_option("-l", "--use-full-matrix", action="store_true",
                      dest="useFullMatrix", help="use full matrix (disable "
                      "removal of dead ends and nonfunctional reactions)")
    parser.set_defaults(tolerance=.95, useFullMatrix=False)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-o")

    if options.tolerance <= 0. or options.tolerance > 1.:
        print("Error: Tolerance must be in interval (0, 1]")
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

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print('\n'+'\n'.join(msgs))

    # 4. Launch flux variability analysis

    print("Launching flux variability analysis.")
    print("tolerance: %g" % options.tolerance)
    if solver == "":
        solver = "default"
    fva = FvAnalyzer(solver)
    minmax, sFlux, lb, ub, ndim = fva.runOnModel(model, fbaParams,
                                                 options.tolerance, options.solutionFile,
                                                 rmDeadEnds=not options.useFullMatrix)

    print ("Info: The reduced network has %u reactions and %u metabolites." %
           ndim[1::-1])
    if not minmax:
        print("Model is infeasible or unbounded. Nothing to do.")
        exit()

    # 5. Write output to file

    try:
        fva.writeSolutionToFile(options.outputFile, model.reactionDict, sFlux,
                                lb, ub, minmax)
    except IOError as strerror:
        print ("Unable to write to file %s:" %
               os.path.basename(options.outputFile))
        print(strerror)
        exit()


if __name__ == "__main__":
    main()
