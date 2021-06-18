#!/usr/bin/env python

""" metano Knockout Analysis

This script performs automatic knockout analysis using Minimization of Metabolic
Adjustment (MOMA) on the metabolic network (given as a set of reaction
equations) and inequality constraints. This is done by sequentially setting the
fluxes through given groups of reactions to zero and performing MOMA against the
wildtype.


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
from builtins import map
from builtins import range
from metano.fba import OptionParser, FbAnalyzer
from metano.moma import MomaAnalyzer, _DEFAULT_ALPHA, _DEFAULT_BETA
from metano.fva import FvAnalyzer
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.defines import (SolverStatus, FbaParam, padNumber, When, decodeWhen,
                     COPYRIGHT_VERSION_STRING)
from numpy import array, dot, nan, inf, empty
import csv, os
try:
    from mpi4py import MPI
    _mpi_avail = True
except ImportError:
    _mpi_avail = False


def _mpi_exit(comm=None, errno=0):
    if comm:
        comm.Abort(errno)
    else:
        exit(errno)

_CUTOFF_LETHAL = 1E-5
_FILE_SUFFIX = '.txt'


def _makeDefaultKoGroups(model):
    """ generate a list of singleton reaction groups from given MetabolicModel
    """
    return [(rea.name, [rea.name]) for rea in model]


def serial_knockout(model, solver, objective, wtSolution, fbaParams,
                    koGroups=None, filePrefix=None, wtObjVal=None, weights=None,
                    numIter=1, useMoma=When.AUTO,
                    lethalityCutoff=_CUTOFF_LETHAL):
    """ perform knockout analysis in serial fashion (single process)

    Keyword arguments:

    model      -- the MetabolicModel
    solver     -- name of QP/NLP solver (or "default")
    objective  -- coefficient vector of linear objective function of FBA
    wtSolution -- FBA solution for wildtype
    fbaParams  -- parameters for FBA (incl. linear equality and inequality
                                      constraints)
    koGroups   -- list of pairs (group name, list of reactions)
    filePrefix -- if given, write MOMA/FBA solutions to files starting with this
                  prefix; if not given, solutions are discarded
    wtObjVal   -- objective function value for wildtype solution (optional)
    weights    -- weight vector for weighted MOMA (None -> perform regular MOMA,
                  else: weight flux i with weights[i])
    numIter    -- number of NLP runs (if NLP solver is used) for MOMA
    useMoma    -- selector of optimization strategy, enum values (class When):
                  NEVER  - always perform FBA (LP)
                  AUTO   - perform FBA first, followed by MOMA if not lethal
                  ALWAYS - always perform MOMA (QP)
    lethalityCutoff
               -- threshold for biomass flux signifying lethality

    Returns list of (distance, diff, obj_val) tuples, indexed like koGroups,

    distance -- value of actual MOMA/FBA objective function at solution
    diff     -- summed absolute difference between mutant and wildtype solutions
    obj_val  -- value of FBA objective function at MOMA/FBA solution for mutant
    """
    if not koGroups:
        koGroups = _makeDefaultKoGroups(model)
#    print "serial: %d reactions, %d knockout groups" % (len(model),
#                                                        len(koGroups))
    wtVec = wtSolution.getVecOrderedByModel(model)
    if wtObjVal is None:
        wtObjVal = dot(objective, wtVec)
    moma = MomaAnalyzer(solver)

    result = []
    for group, reaList in koGroups:
        print(group)

        # Restrict flux through all reactions in group to zero
        tmp_lb = {}  # {index in model : LB value}
        tmp_ub = {}  # {  - '' -       : UB  ''  }
        for rea in reaList:
            try:
                rIndex = model.reactionDict[rea]
            except KeyError:
                # Skip any blocked reactions (knockout has no effect)
                continue
            tmp_lb[rIndex] = model.reactions[rIndex].lb
            tmp_ub[rIndex] = model.reactions[rIndex].ub
            model.reactions[rIndex].lb = model.reactions[rIndex].ub = 0.

        if not tmp_lb:
            print("  -> skipped (only blocked reactions)")
            result.append((0., 0., wtObjVal))
            if filePrefix:
                wtSolution.writeToFile(filePrefix+group+_FILE_SUFFIX)
            continue

        obj_val = lethalityCutoff + 1.

        # Perform FBA first to check for lethality
        if useMoma != When.ALWAYS:
            fba = FbAnalyzer(fbaParams.solver)
            # Note: Model is already reduced
            obj_val, solutionFlux = fba.runOnModel(model, fbaParams,
                                                   rmDeadEnds=False)[:2]
            status = SolverStatus.PREPARED
            dist = None
            if not solutionFlux:
                obj_val = 0.

        # Perform MOMA if knockout not already predicted by FBA to be lethal
        if useMoma != When.NEVER and abs(obj_val) > lethalityCutoff:
            dist, solutionFlux, status = moma.runOnModel(model, wtSolution,
                fbaParams.linConstraints, numIter, weights)[:3]

        if status == SolverStatus.UNKNOWN:
            result.append((nan, nan, 0.))
        elif status in (SolverStatus.PRIM_INFEAS, SolverStatus.DUAL_INFEAS):
            result.append((inf, inf, 0.))
        elif len(solutionFlux) == 0:
            if status == SolverStatus.PREPARED:
                # FBA infeasible
                result.append((inf, inf, 0.))
            else:
                result.append((nan, nan, 0.))
        else:
            diff = sum(solutionFlux.absDiff(wtSolution).fluxDict.values())
            if dist is None:
                dist = moma.evalObjFunc(wtVec,
                    solutionFlux.getVecOrderedByModel(model), weights)
            else:
                obj_val = max(0., dot(objective,
                                      solutionFlux.getVecOrderedByModel(model)))
            result.append((dist, diff, obj_val))
            if filePrefix:
                solutionFlux.writeToFile(filePrefix+group+_FILE_SUFFIX)

        # Restore original constraints
        for rea in reaList:
            try:
                rIndex = model.reactionDict[rea]
            except KeyError:
                continue
            model.reactions[rIndex].lb = tmp_lb[rIndex]
            model.reactions[rIndex].ub = tmp_ub[rIndex]

    return result


def parallel_ko_dispatch(comm, numprocs, model, solver, objective, wtSolution,
                         wtObjVal=None, koGroups=None, filePrefix=None,
                         weights=None):
    """ root process (dispatcher) for parallel knockout analysis

    Keyword arguments:

    comm       -- MPI communicator
    numprocs   -- number of processes (including root process)
    model      -- the MetabolicModel
    solver     -- name of solver (or "default")
    objective  -- coefficient vector of linear objective function of FBA
    wtSolution -- FBA solution for wildtype
    wtObjVal   -- objective function value for wildtype solution (optional)
    koGroups   -- list of pairs (group name, list of reactions)
    filePrefix -- if given, write MOMA/FBA solutions to files starting with this
                  prefix; if empty, solutions are discarded
    weights    -- weight vector for weighted MOMA (None -> perform regular MOMA,
                  else: weight flux i with weights[i])

    Returns list of (distance, diff, obj_val) tuples, indexed like koGroups,

    distance -- value of actual MOMA/FBA objective function at solution
    diff     -- summed absolute difference between mutant and wildtype solutions
    obj_val  -- value of FBA objective function at MOMA/FBA solution for mutant
    """
    if not koGroups:
        koGroups = _makeDefaultKoGroups(model)
#    print "dispatch: %d reactions, %d knockout groups" % (len(model),
#                                                          len(koGroups))
    wtVec = wtSolution.getVecOrderedByModel(model)
    if wtObjVal is None:
        wtObjVal = dot(objective, wtVec)
    moma = MomaAnalyzer(solver)

    nJobs = len(koGroups)
    stopMsg = array(-1, 'i')  # stop message (-1)
    status = MPI.Status()
    result = [None]*nJobs
    jobList = [-1]*numprocs  # mapping of process -> job (jobList[0] is dummy)
    # receive buffer for MOMA/FBA solution vectors
    solutionVec = empty(len(model), 'd')
    # index of next job to be processed (as 0-dimensional numpy array)
    nxt = array(0, 'i')
    free_processes = set(range(1, numprocs))  # set of processes free for work
    nDone = 0  # number of finished jobs

    while nDone < nJobs:
        while free_processes and nxt < nJobs:
            # Assign reaction to free process
            print(koGroups[nxt][0])
            proc_index = free_processes.pop()
            # Send index of next job to next free process
            comm.Send(nxt, proc_index)
            # Mark in job list (int() is used to convert from ref. to literal)
            jobList[proc_index] = int(nxt)
            nxt += 1

        # If there are more processes than jobs, stop extra processes
        while free_processes:
            comm.Send(stopMsg, free_processes.pop())

        # Receive results from compute processes
        comm.Recv(solutionVec, MPI.ANY_SOURCE, status=status)
        proc_index = status.Get_source()

        # Test if solutionVec[0] is nan (nan means not converged)
        if solutionVec[0] != solutionVec[0]:
            result[jobList[proc_index]] = (nan, nan, 0.)
        elif solutionVec[0] == inf:
            # inf means infeasible, i.e. no solution exists
            result[jobList[proc_index]] = (inf, inf, 0.)
        else:
            jobIndex = jobList[proc_index]
            solutionFlux = MetabolicFlux(model, solutionVec)
            diff = sum(solutionFlux.absDiff(wtSolution).fluxDict.values())
            obj_val = max(0., dot(objective, solutionVec))
            distance = moma.evalObjFunc(wtVec, solutionVec, weights)
            result[jobIndex] = (distance, diff, obj_val)
            if filePrefix:
                solutionFlux.writeToFile(filePrefix+koGroups[jobIndex][0]+
                                         _FILE_SUFFIX)

        # Increase job counter and mark reporting process as free
        nDone += 1
        free_processes.add(proc_index)

    # Stop last compute process
    while free_processes:
        comm.Send(stopMsg, free_processes.pop())

    return result


def parallel_ko_worker(comm, model, solver, wtVec, fbaParams, koGroups=None,
                       weights=None, numIter=1, useMoma=When.AUTO,
                       lethalityCutoff=_CUTOFF_LETHAL):
    """ compute process (worker) for parallel knockout analysis

    Keyword arguments:

    comm       -- MPI communicator
    model      -- the MetabolicModel
    solver     -- name of QP/NLP solver (or "default")
    wtVec      -- FBA solution for wildtype
    fbaParams  -- parameters for FBA (incl. linear equality and inequality
                                      constraints)
    koGroups   -- list of pairs (group name, list of reactions)
    weights    -- weight vector for weighted MOMA (None -> perform regular MOMA,
                  else: weight flux i with weights[i])
    numIter    -- number of NLP runs (if NLP solver is used) for MOMA
    useMoma    -- selector of optimization strategy, enum values (class When):
                  NEVER  - always perform FBA (LP)
                  AUTO   - perform FBA first, followed by MOMA if not lethal
                  ALWAYS - always perform MOMA (QP)
    lethalityCutoff
               -- threshold for biomass flux signifying lethality
    """
    if not koGroups:
        koGroups = _makeDefaultKoGroups(model)
#    print "worker: %d reactions, %d knockout groups" % (len(model),
#                                                        len(koGroups))
    wtSolution = MetabolicFlux(model, wtVec)
    nanVec = array([nan]*len(wtVec))
    infVec = array([inf]*len(wtVec))
    moma = MomaAnalyzer(solver)

    # index of reaction group to be knocked out (as 0-dimensional numpy array)
    i = array(0, 'i')
    comm.Recv(i)
    while i >= 0:
        group, reaList = koGroups[i]

        # Restrict flux through all reactions in group to zero
        tmp_lb = {}  # {index in model : LB value}
        tmp_ub = {}  # {  - '' -       : UB  ''  }
        for rea in reaList:
            try:
                rIndex = model.reactionDict[rea]
            except KeyError:
                # Skip any blocked reactions (knockout has no effect)
                continue
            tmp_lb[rIndex] = model.reactions[rIndex].lb
            tmp_ub[rIndex] = model.reactions[rIndex].ub
            model.reactions[rIndex].lb = model.reactions[rIndex].ub = 0.

        if not tmp_lb:
            print("  - skipped group '%s' (only blocked reactions)" % group)
            comm.Send(wtVec)
            comm.Recv(i)
            continue

        obj_val = lethalityCutoff + 1.

        # Perform FBA first to check for lethality
        if useMoma != When.ALWAYS:
            fba = FbAnalyzer(fbaParams.solver)
            # Note: Model is already reduced
            obj_val, solutionFlux = fba.runOnModel(model, fbaParams,
                                                   rmDeadEnds=False)[:2]
            status = SolverStatus.PREPARED
            if not solutionFlux:
                obj_val = 0.

        # Perform MOMA if knockout not already predicted by FBA to be lethal
        if useMoma != When.NEVER and abs(obj_val) > lethalityCutoff:
            solutionFlux, status = moma.runOnModel(model, wtSolution,
                fbaParams.linConstraints, numIter, weights)[1:3]

        # Send result to dispatcher
        if status == SolverStatus.UNKNOWN:
            # If status is 'unknown, not converged', send vector with only 'nan'
            # entries
            comm.Send(nanVec)
        elif status in (SolverStatus.PRIM_INFEAS, SolverStatus.DUAL_INFEAS):
            # If status is 'infeasible', send vector with only 'inf' entries
            comm.Send(infVec)
        elif len(solutionFlux) == 0:
            if status == SolverStatus.PREPARED:
                # FBA infeasible
                comm.Send(infVec)
            else:
                comm.Send(nanVec)
        else:
            solutionVec = array(solutionFlux.getVecOrderedByModel(model))
            comm.Send(solutionVec)

        # Restore original constraints
        for rea in reaList:
            try:
                rIndex = model.reactionDict[rea]
            except KeyError:
                continue
            model.reactions[rIndex].lb = tmp_lb[rIndex]
            model.reactions[rIndex].ub = tmp_ub[rIndex]

        # Get next index
        comm.Recv(i)


def main():
    if _mpi_avail:
        comm = MPI.COMM_WORLD
        numprocs = comm.Get_size()
        rank = comm.Get_rank()
        if numprocs == 1:
            comm = None
    else:
        comm = None

    if not _mpi_avail or rank == 0:    # Start of root process code

        # 1. Parse command line

        usage = "Usage: %prog [options]"
        version = "Knockout analysis\n" + COPYRIGHT_VERSION_STRING
        parser = OptionParser(usage=usage, version=version)
        momaOptions = parser.add_option_group("General MOMA options")
        momaOptions.add_option("-r", "--reactions", dest="reactionFile",
                               help="perform knockout analysis on the network "
                               "given by the reaction FILE", metavar="FILE")
        momaOptions.add_option("-p", "--parameters", dest="paramFile",
                               help="use the given scenario FILE for the "
                               "wildtype", metavar="FILE")
        momaOptions.add_option("-o", "--output", dest="outputFile", help="write"
                               " output of knockout analysis to FILE",
                               metavar="FILE")
        momaOptions.add_option("-w", "--wt-solution", dest="wtSolution",
                               help="use the given solution FILE for the "
                               "wildtype (optional)", metavar="FILE")
        momaOptions.add_option("-v", "--reduce-by-fva", dest="redFvaFile",
                               help="use FVA result from FILE to reduce matrix "
                               "(experimental; only works if perturbed solution"
                               " space is a subset of the wildtype solution "
                               "space)", metavar="FILE")
        koOptions = parser.add_option_group("Options for knockout analysis")
        koOptions.add_option("-f", "--solution-files", dest="fluxFilePrefix",
                             help="write MOMA/FBA solutions to files marked "
                             "with FILE-PREFIX (may contain path)",
                             metavar="FILE-PREFIX")
        koOptions.add_option("-g", "--ko-groups", dest="koGroupDictFile",
                             help="read reaction groups to be knocked out "
                             "together from CSV FILE (if not given, individual "
                             "reactions are knocked out)", metavar="FILE")
        koOptions.add_option("-n", "--ko-ungrouped", action="store_true",
                             dest="koUngrouped", help="knock out individual "
                             "reactions that do not belong to any group (only "
                             "with -g)")
        koOptions.add_option("-m", "--moma", dest="useMoma", help="perform MOMA"
                             ", options: 'always', 'never' (FBA only), 'auto' "
                             "(default: FBA first, then MOMA if not lethal)",
                             metavar="WHEN")
        wmOptions = parser.add_option_group("Extra options for weighted MOMA")
        wmOptions.add_option("-x", "--wmoma", dest="wMomaFvaFile",
                             help="perform weighted MOMA with weights computed "
                             "from the given FVA solution FILE", metavar="FILE")
        wmOptions.add_option("-a", "--alpha", dest="alpha", type="float",
                             help="parameter ALPHA for computation of weights:"
                             " w = alpha + exp(-beta*(fvaMax-fvaMin)) (default:"
                              " %g)" % _DEFAULT_ALPHA)
        wmOptions.add_option("-b", "--beta", dest="beta", type="float",
                             help="parameter BETA for computation of weights: "
                             "w = alpha + exp(-beta*(fvaMax-fvaMin)) (default: "
                             "%g)" % _DEFAULT_BETA)
        optOptions = parser.add_option_group("Options for optimization")
        optOptions.add_option("-s", "--solver", dest="solver", help="QP solver "
                              "to be used for MOMA (default cvxopt)",
                              metavar="SOLVER")
        optOptions.add_option("-i", "--iterations", dest="numIter", type="int",
                              help="number of NLP runs to perform (from "
                              "different random start points) - per knockout")
        optOptions.add_option("-l", "--use-full-matrix", action="store_true",
                              dest="useFullMatrix", help="use full matrix "
                              "(disable removal of dead ends and nonfunctional "
                              "reactions)")
        parser.set_defaults(solver="default", numIter=1, koUngrouped=False,
                            useMoma="auto", useFullMatrix=False,
                            alpha=_DEFAULT_ALPHA, beta=_DEFAULT_BETA)

        try:
            options, _ = parser.parse_args()
            parser.check_required("-r")
            parser.check_required("-p")
            parser.check_required("-o")
        except SystemExit as e:
            _mpi_exit(comm, e.args[0])

        weighted = bool(options.wMomaFvaFile)
        if weighted:
            if options.alpha < 0.:
                print("Error: alpha must be non-negative.")
                _mpi_exit(comm, 1)
            if options.beta <= 0.:
                print("Error: beta must be positive.")
                _mpi_exit(comm, 1)

        if options.numIter < 1:
            print("Error: Number of NLP runs must be positive.")
            _mpi_exit(comm, 1)

        try:
            useMoma = decodeWhen(options.useMoma, set((When.NEVER, When.AUTO,
                                                       When.ALWAYS)))
        except KeyError:
            print ("Error: Illegal value for -m option (allowed values are "
                   "'auto', 'never', or 'always').")
            _mpi_exit(comm, 1)

        # 2. Parse reaction file

        rparser = ReactionParser()
        model = MetabolicModel()
        try:
            model.addReactionsFromFile(options.reactionFile, rparser)
        except IOError as strerror:
            print ("An error occurred while trying to read file %s:" %
                   os.path.basename(options.reactionFile))
            print(strerror)
            _mpi_exit(comm, 1)
        except SyntaxError as strerror:
            print ("Error in reaction file %s:" %
                   os.path.basename(options.reactionFile))
            print(strerror)
            _mpi_exit(comm, 1)

        # 3. Parse scenario file

        model_messages = []
        pparser = ParamParser()
        try:
            maxmin, objStr, fbaSolver, numIter, lb , ub = \
                pparser.parse(options.paramFile)
            # Set flux bounds in model
            model.setFiniteBounds(lb, ub, True, model_messages)
            if fbaSolver == "":
                fbaSolver = "default"
            fbaParams = FbaParam(fbaSolver, maxmin, objStr, numIter)
            lin_constraints = pparser.lin_constraints
            fbaParams.setLinConstraints(lin_constraints)
        except IOError as strerror:
            print ("An error occurred while trying to read file %s:" %
                   os.path.basename(options.paramFile))
            print(strerror)
            _mpi_exit(comm, 1)
        except SyntaxError as strerror:
            print ("Error in scenario file %s:" %
                   os.path.basename(options.paramFile))
            print(strerror)
            _mpi_exit(comm, 1)
        except ValueError as strerror:
            print(strerror)
            _mpi_exit(comm, 1)

        # Show warning and info messages of parsers
        msgs = (rparser.getMessages() + pparser.getMessages() +
                [x[1] for x in model_messages])
        if msgs:
            print('\n'+'\n'.join(msgs))

        # 4. Read groups of reactions from file (if given)

        koGroups = []  # list of pairs (group name, list of reactions)

        if options.koGroupDictFile:
            koGroupDict = {}     # dict {group name : list of reactions}
            koGroupDictRev = {}  # dict {reaction : list of groups}
            koOrder = []         # order of reactions
            try:
                with open(options.koGroupDictFile) as f:
                    csvReader = csv.reader(f, delimiter=';')

                    for row in csvReader:
                        group, rea = list(map(str.strip, row[:2]))
                        koOrder.append(rea)

                        if group in koGroupDict:
                            koGroupDict[group].append(rea)
                        else:
                            koGroupDict[group] = [rea]

                        if rea in koGroupDictRev:
                            koGroupDictRev[rea].append(group)
                        else:
                            koGroupDictRev[rea] = [group]

            except csv.Error as e:
                print ("Error in CSV file %s, line %d: %s" %
                       (os.path.basename(options.koGroupDictFile),
                        csvReader.line_num, e))
                _mpi_exit(comm, 1)
            except ValueError as strerror:
                print ("Error in CSV file %s, line %d: %s" %
                       (os.path.basename(options.koGroupDictFile),
                        csvReader.line_num, strerror))
                _mpi_exit(comm, 1)
            except IOError as strerror:
                print ("An error occurred while trying to read file %s:" %
                       os.path.basename(options.koGroupDictFile))
                print(strerror)
                _mpi_exit(comm, 1)

            # Generate list of knockouts (order as close as possible to order of
            #                             reactions in file)
            processedGroups = set()
            for rea in koOrder:
                for group in koGroupDictRev[rea]:

                    if group in ("", "NULL"):
                        if options.koUngrouped:
                            # Create new group "KO_x" for individual reaction x
                            koGroups.append(("KO_"+rea, [rea]))
                        # else skip

                    elif group not in processedGroups:
                        koGroups.append((group, koGroupDict[group]))
                        processedGroups.add(group)

        else:
            # Knock out individual reactions
            koGroups = _makeDefaultKoGroups(model)

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
                _mpi_exit(comm, 1)
            except SyntaxError as strerror:
                print ("Error in FVA solution file %s:" %
                       os.path.basename(options.redFvaFile))
                print(strerror)
                _mpi_exit(comm, 1)

            if not fbaSolution.hasSameReactions(model):
                print ("Error in file %s: FVA solution and model must have the "
                       "same reactions." % os.path.basename(options.redFvaFile))
                _mpi_exit(comm, 1)

            # Determine blocked reactions based on FVA solution
            blockedReactions = FvAnalyzer.getBlockedReactions(redMinmax)

        elif options.useFullMatrix:
            blockedReactions = []
        else:
            # Determine blocked reactions based on dead-end analysis
            blockedReactions = model.findDeadEnds(True)[1]

        # 6. (Optionally) Read FVA solution and compute weights for wMOMA

        if weighted:
            try:
                fbaSolution, wmomaMinmax = \
                    FvAnalyzer.parseSolutionFile(options.wMomaFvaFile)
            except IOError as strerror:
                print ("An error occurred while trying to read file %s:" %
                       os.path.basename(options.wMomaFvaFile))
                print(strerror)
                _mpi_exit(comm, 1)
            except SyntaxError as strerror:
                print ("Error in FVA solution file %s:" %
                       os.path.basename(options.wMomaFvaFile))
                print(strerror)
                _mpi_exit(comm, 1)

            if not fbaSolution.hasSameReactions(model):
                print ("Error in file %s: FVA solution and model must have the "
                     "same reactions." % os.path.basename(options.wMomaFvaFile))
                _mpi_exit(comm, 1)

            weights = MomaAnalyzer.getWeightsFromFluxVar(model, wmomaMinmax,
                                                    options.alpha, options.beta)

        # 7.a Read wildtype solution from file (if given)

        if options.wtSolution:
            # Read wildtype solution from file
            wtSolution = MetabolicFlux()
            try:
                wtSolution.readFromFile(options.wtSolution)
            except IOError as strerror:
                print ("An error occurred while trying to read file %s:" %
                       os.path.basename(options.wtSolution))
                print(strerror)
                _mpi_exit(comm, 1)
            except SyntaxError as strerror:
                print ("An error occurred parsing file %s:" %
                       os.path.basename(options.wtSolution))
                print(strerror)
                _mpi_exit(comm, 1)

            if not wtSolution.hasSameReactions(model):
                print ("Error in file %s: FBA solution and model must have the "
                       "same reactions." % os.path.basename(options.wtSolution))
                _mpi_exit(comm, 1)

            # Build coefficient vector of linear objective function
            objective = array(ParamParser.convertObjFuncToLinVec(objStr,
                                                            model.reactionDict))
            wt_obj_value = dot(objective,
                               wtSolution.getVecOrderedByModel(model))

        # 7.b If FBA solution is not given, perform FBA

        else:
            # Compute wildtype solution via FBA
            fba = FbAnalyzer(fbaSolver)

            wt_obj_value, wtSolution, ndim = fba.runOnModel(model, fbaParams,
                rmDeadEnds=not options.useFullMatrix)
            if not options.useFullMatrix:
                print ("Info: The reduced network for FBA has %u reactions and "
                       "%u metabolites." % ndim[1::-1])
            if len(wtSolution) == 0:
                print("Model is infeasible or unbounded. Nothing to do.")
                _mpi_exit(comm)

        solver = options.solver
        if not solver:
            solver = "default"
        numIter = options.numIter

        # 8. Generate submodel with blocked reactions removed

        if blockedReactions:
            modelRed = model.getSubModelByExcludeList(blockedReactions)
            reactionsRed = set(modelRed.getReactionNames())
            # Reduce solution to non-dead reactions
            wtSolutionRed = MetabolicFlux()
            for rea in wtSolution:
                if rea in reactionsRed:
                    wtSolutionRed.fluxDict[rea] = wtSolution.fluxDict[rea]
                    wtSolutionRed.boundsDict[rea] = wtSolution.boundsDict[rea]
            # Also reduce weights vector
            if weighted:
                weightsRed = [0.]*len(modelRed)
                for rea in wtSolution:
                    if rea in reactionsRed:
                        weightsRed[modelRed.reactionDict[rea]] = \
                            weights[model.reactionDict[rea]]
                weightVec = array(weightsRed)
            else:
                weightVec = None
        else:
            modelRed = model
            wtSolutionRed = wtSolution
            if weighted:
                weightVec = array(weights)
            else:
                weightVec = None

        modelStr = modelRed.writeToString(False, False, False)
        nReactions = len(modelRed)
        lb, ub = list(map(array, modelRed.getBounds()))
        # Build coefficient vector of linear objective function
        objective = array(ParamParser.convertObjFuncToLinVec(objStr,
                                                        modelRed.reactionDict))
        if _mpi_avail and numprocs > 1:
            intBuf = array((nReactions, numIter, int(useMoma), int(maxmin),
                            int(weighted)), 'i')
            comm.Bcast(intBuf)
            wtVec = array(wtSolutionRed.getVecOrderedByModel(modelRed))
            if weighted:
                weightVec = array(weightsRed)
            else:
                weightVec = None

    else:                              # Start of compute process (slave) code
        solver = None
        modelStr = None
        intBuf = empty(5, 'i')
        comm.Bcast(intBuf)
        nReactions, numIter, useMoma, maxmin, weighted = intBuf
        lb = empty(nReactions, 'd')
        ub = empty(nReactions, 'd')
        wtVec = empty(nReactions, 'd')
        fbaSolver = None
        objStr = None
        lin_constraints = None
        koGroups = None
        if weighted:
            weightVec = empty(nReactions, 'd')
        else:
            weightVec = None

    # 9. Sequentially knock out fluxes and perform MOMA against the wildtype
    #    solution

    if _mpi_avail and numprocs > 1:
        # Launch parallel analysis
        solver = comm.bcast(solver)
        modelStr = comm.bcast(modelStr)
        comm.Bcast(lb)
        comm.Bcast(ub)
        comm.Bcast(wtVec)
        fbaSolver = comm.bcast(fbaSolver)
        objStr = comm.bcast(objStr)
        fbaParams = FbaParam(fbaSolver, maxmin, objStr)
        fbaParams.linConstraints = comm.bcast(lin_constraints)
        koGroups = comm.bcast(koGroups)
        if weighted:
            comm.Bcast(weightVec)

        if rank == 0:
            result = parallel_ko_dispatch(comm, numprocs, modelRed, solver,
                            objective, wtSolutionRed, wt_obj_value, koGroups,
                            options.fluxFilePrefix, weightVec)
        else:
            model = MetabolicModel()
            model.addReactionsFromString(modelStr)
            model.setBounds(lb, ub)
            parallel_ko_worker(comm, model, solver, wtVec, fbaParams, koGroups,
                               weightVec, numIter, useMoma)

    else:
        # Launch serial analysis
        result = serial_knockout(modelRed, solver, objective, wtSolutionRed,
                                 fbaParams, koGroups, options.fluxFilePrefix,
                                 wt_obj_value, weightVec, numIter, useMoma)

    if not _mpi_avail or rank == 0:

    # 10. Write output to file

        print("Writing output file.")
        try:
            write_output(options.outputFile, [x[0] for x in koGroups], result,
                         wt_obj_value)
        except IOError as strerror:
            print ("Unable to write to file %s:" %
                   os.path.basename(options.outputFile))
            print(strerror)
            _mpi_exit(comm, 1)


def write_output(filename, group_names, result, wt_obj_val):
    # Construct output table and get widths of table columns (length of
    # longest occurring string)
    output = []
    for i in range(len(group_names)):
        output.append((group_names[i],) + tuple(map(padNumber, list(map(repr,
            (100.*result[i][2]/wt_obj_val,) + result[i])))))
        maxlen = tuple(map(len, output[-1]))

    # Create list of (index, value of objective function) pairs, sorted by value
    # in descending order - for ordering the rows of the output table
    order = sorted(zip(list(range(len(result))), [x[2] for x in result]),
                   key=lambda x : x[1], reverse=True)

    with open(filename, 'w') as f:
        if filename.split(".")[-1] == "csv":
            csvwriter = csv.writer(f, delimiter=";", quotechar='"', quoting=csv.QUOTE_MINIMAL)
            # Write Header to CSV file
            csvwriter.writerow(["NAME", "OBJ_FUN(%)", "DISTANCE", "ABS_DIFF", "OBJ_FUN_VAL"])
            # Construct output table and get widths of table columns (length of
            # longest occurring string)
            for i, _ in order:
                reaName, objfPercent, dist, diff, objfValue = output[i]
                csvwriter.writerow(output[i])

        else:
            # Write table head
            f.write("NAME".ljust(maxlen[0])+"   "+"OBJ_FUN(%)".center(maxlen[1])+" "
                    +"DISTANCE".rjust(maxlen[2])+" "+"ABS_DIFF".rjust(maxlen[3])+" "
                    +"OBJ_FUN_VAL".rjust(maxlen[4])+"\n")

            # Write table rows
            for i, _ in order:
                reaName, objfPercent, dist, diff, objfValue = output[i]
                f.write(reaName.ljust(maxlen[0]) + " : " +
                        objfPercent.rjust(maxlen[1]) + " " + dist.rjust(maxlen[2]) +
                        " " + diff.rjust(maxlen[3]) + " " +
                        objfValue.rjust(maxlen[4]) + "\n")


if __name__ == "__main__":
    main()
