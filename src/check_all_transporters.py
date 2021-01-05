#!/usr/bin/env python

""" metano Script for checking transporters

This script checks all meaningful combinations of importers (e.g. C source,
N source, S source). It supports auxiliary importers of the same "type", e.g.
acetone import + glucose import as auxiliary transporter.


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
from metano.defines import FbaParam, COPYRIGHT_VERSION_STRING
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.transporterparser import TransporterParser
from metano.fba import OptionParser, FbAnalyzer
from numpy import array, inf
import os

# debug flag: if True, the LB/UB values for every FBA are dumped to the console
debugFlag = False


def checkRecursive(resultList, currentList, otherLists, fba, matrix, reactions,
                   reaction_names, lb, ub, fbaParams, partialResult=[]):
    """ recursively enumerate all meaningful combinations of transporters and
        launch FBA on each such combination

    Keyword arguments:

    resultList     -- list for storage of results on innermost recursion level
    currentList    -- list of importers of the same type (perform FBA if empty)
    otherLists     -- list of lists of importers of other types (may be empty)
    fba            -- FbAnalyzer to be used for FBA
    matrix         -- the stoichiometric matrix of the metabolic network
    reactions      -- dictionary of all reactions { name : matrix column }
    reaction_names -- list of reaction names (column index -> name)
    lb             -- list of lower bounds (indexed like matrix columns)
    ub             -- list of upper bounds (indexed like matrix columns)
    fbaParams      -- ParamParser object for getting nonlinear constraints
    partialResult  -- growing result dictionary

    Returns: nothing, appends results to resultList

    resultList is list of dictionaries containing keys <type>, <type>_co,
    <type>_flux, <type>_co_flux for every importer and co-importer and key
    "solution" for the FBA solution
    """
    if currentList == []:
        # Recursion anchor: We have picked transporters of every type -> run FBA

        matrixSplit, reactionsSplit, lbSplit, ubSplit = \
            FbAnalyzer.splitFluxes(matrix, reaction_names, lb, ub)

        solutionSplit = fba.run(reactionsSplit, matrixSplit, lbSplit, ubSplit,
                                fbaParams)[1]

        keys = list(partialResult.keys())
        if solutionSplit == []:
            solution = []
            for key in keys:
                rea = partialResult[key]
                partialResult[key+"_flux"] = None
        else:
            solution = array(FbAnalyzer.rejoinFluxes(solutionSplit,
                                                     reactionsSplit,
                                                     reactions))
            for key in keys:
                rea = partialResult[key]
                if rea in reactions:
                    partialResult[key+"_flux"] = solution[reactions[rea]]
                else:
                    partialResult[key+"_flux"] = None

        partialResult["solution"] = solution

        # for debugging: store LB/UB as dictionary
        if debugFlag:
            lbub = {}
            for rea in reactions:
                rea_index = reactions[rea]
                lbub[rea] = (lb[rea_index], ub[rea_index])
            partialResult["bounds"] = lbub

        resultList.append(partialResult)
    else:
        # It's sufficient to query the first element as all members of
        # currentList[i] have the same type.
        typ = currentList[0][0]

        # Set next currentList for recursion (head of otherLists if not empty)
        if otherLists == []:
            nextList = []
        else:
            nextList = otherLists[0]

        for (_, trans_name, factor, cotrans) in currentList:
            # Make a copy of partialResult and extend by keys for current type
            tmpResult = dict(partialResult)
            tmpResult[typ] = trans_name
            tmpResult[typ+"_co"] = None

            if factor is not None:
                # Set flux through transporter to fix value
                index = reactions[trans_name]
                tmp_lb_trans = lb[index]
                tmp_ub_trans = ub[index]
                lb[index] = ub[index] = factor
                if cotrans is not None:
                    # make another copy of partial result
                    tmpResBak = dict(tmpResult)

                # First perform FBA without cotransporter
                checkRecursive(resultList, nextList, otherLists[1:], fba,
                               matrix, reactions, reaction_names, lb, ub,
                               fbaParams, tmpResult)

                # Now set flux through cotransporter to appropriate value
                # (factor)
                if cotrans is not None:
                    # get copy of partial result
                    tmpResult = tmpResBak
                    tmpResult[typ+"_co"] = cotrans

                    # Find factor cotrans_fac in currentList
                    for (_, cotrans_name, cotrans_fac, _) in currentList:
                        if cotrans_name == cotrans:
                            break

                    co_index = reactions[cotrans]
                    tmp_lb_cotrans = lb[co_index]
                    tmp_ub_cotrans = ub[co_index]
                    lb[co_index] = ub[co_index] = cotrans_fac

                    checkRecursive(resultList, nextList, otherLists[1:], fba,
                                   matrix, reactions, reaction_names, lb, ub,
                                   fbaParams, tmpResult)

                    lb[co_index] = tmp_lb_cotrans
                    ub[co_index] = tmp_ub_cotrans

                lb[index] = tmp_lb_trans
                ub[index] = tmp_ub_trans
            else:
                # Allow arbitrary flux through transporter
                index = reactions[trans_name]
                tmp_lb_trans = lb[index]
                tmp_ub_trans = ub[index]
                lb[index] = -inf
                ub[index] = inf

                # First perform FBA without cotransporter
                checkRecursive(resultList, nextList, otherLists[1:], fba,
                               matrix, reactions, reaction_names, lb, ub,
                               fbaParams, tmpResult)

                # Now also allow arbitrary flux through cotransporter
                if cotrans is not None:
                    # make another copy of partial result
                    tmpResult = dict(tmpResult)
                    tmpResult[typ+"_co"] = cotrans
                    co_index = reactions[cotrans]
                    tmp_lb_cotrans = lb[co_index]
                    tmp_ub_cotrans = ub[co_index]
                    lb[co_index] = -inf
                    ub[co_index] = inf

                    checkRecursive(resultList, nextList, otherLists[1:], fba,
                                   matrix, reactions, reaction_names, lb, ub,
                                   fbaParams, tmpResult)

                    lb[co_index] = tmp_lb_cotrans
                    ub[co_index] = tmp_ub_cotrans

                lb[index] = tmp_lb_trans
                ub[index] = tmp_ub_trans


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "Transporter check\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform Flux Balance Analysis on the network given "
                           "by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use the given scenario FILE for Flux Balance "
                           "Analysis", metavar="FILE")
    parser.add_option("-t", "--transporters", dest="transFile",
                      help="use list of transporters from FILE", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write output of Flux Balance Analysis to FILE",
                      metavar="FILE")

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-t")
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
    reactions = model.reactionDict
    matrix = array(model.getStoichiometricMatrix())

    # 3. Parse scenario file

    model_messages = []
    pparser = ParamParser()
    try:
        # Parse file, get maxmin, name of objective function, and solver name
        maxmin, obj_name, solver, numIter, lb, ub = pparser.parse(
            options.paramFile)
        fbaParams = FbaParam(solver, maxmin, obj_name, numIter)
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
    lb, ub = list(map(array, model.getBounds()))

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print('\n'+'\n'.join(msgs))

    # 4. Parse transporter file

    tparser = TransporterParser()
    try:
        transporters = tparser.parse(options.transFile)
        transpByType = tparser.split_by_type()
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.transFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print(strerror)
        exit()
    except ValueError as strerror:
        print(strerror)
        exit()

    # 5. Set all transporter fluxes to zero

    for (trans_name, _, _) in transporters:
        if trans_name not in reactions:
            print ("Error: Transporter %s is not a reaction of the network." %
                   trans_name)
            exit()
        index = reactions[trans_name]
        lb[index] = ub[index] = 0.

    # 6. Identify main ("carbon") source - the one that has factors

    main_type = None
    for entry in transpByType:
        typ = entry[0][0]
        index = 0
        for (_, _, factor, _) in entry:

            if main_type is None and factor is not None:
                if index != 0:
                    print ("Error: Main transporter type must have flux "
                           "factors. (Identified main type '%s')" % main_type)
                    exit()
                main_type = typ

            elif typ == main_type and factor is None:
                print ("Error: Main transporter type must have flux "
                       "factors. (Identified main type '%s')" % main_type)
                exit()

            elif typ != main_type and factor is not None:
                print ("Error: Only one transporter type can have factors.\n"
                       "Found factors with types %s and %s" % (main_type, typ))
                exit()

            index += 1

    for i in range(len(transpByType)):
        if transpByType[i][0][0] == main_type:
            main_index = i
            break

    # 7. Perform flux balance analysis for all combinations of transporters

    if solver == "":
        solver = "default"
    fba = FbAnalyzer(solver)

    resultList = []
    checkRecursive(resultList, transpByType[main_index],
                   transpByType[:main_index] +
                   transpByType[main_index+1:], fba, matrix,
                   reactions, model.getReactionNames(), lb, ub, fbaParams)

    # For debugging: Output lb/ub to console
    if debugFlag:
        keys = list(resultList[0].keys())
        keys.remove("solution")
        keys.remove("bounds")
        keys.sort()
        for result in resultList:
            print("***", end=' ')
            for key in keys:
                print(key+":", result[key], end=' ')
            print()
            lbub = result["bounds"]
            for rea in sorted(lbub):
                print(rea, lbub[rea])

    # 8. Write output to file

    try:
        write_output(options.outputFile, reactions, resultList, main_type)
    except IOError as strerror:
        print ("Unable to write to file %s:" %
               os.path.basename(options.outputFile))
        print(strerror)
        exit()


def write_output(filename, reactions, resultList, main_type):
    if resultList == []:
        return  # Nothing to write

    # Get all types
    types = []
    for key in resultList[0]:
        try:
            if not key.endswith("_co") and resultList[0][key] in reactions:
                types.append(key)
        except TypeError:
            # Skip non-hashable values of resultList[0][key], such as dicts or
            # lists
            pass

    # Get length of longest reaction name (for formatting)
    maxtypelen = len(max(types, key=len))+len("_cosource")
    maxlen = max(len(max(reactions, key=len)), maxtypelen)

    # Determine optimal column widths
    col_width = []
    for i in range(len(resultList)):
        width = 12
        for typ in types:
            width = max(width, len(resultList[i][typ]))
            if resultList[i][typ+"_co"] in reactions:
                width = max(width, len(resultList[i][main_type+"_co"]))
        col_width.append(width)

    with open(filename, 'w') as f:
        # First line: transporter - main_type first
        line = (main_type+"_source").ljust(maxlen)+" :"
        for i in range(len(resultList)):
            line += ' ' + resultList[i][main_type].rjust(col_width[i])
        f.write(line+'\n')

        # Second line: cotransporter
        line = (main_type+"_cosource").ljust(maxlen)+" :"
        any_cotrans = False
        for i in range(len(resultList)):
            if resultList[i][main_type+"_co"] in reactions:
                s = resultList[i][main_type+"_co"]
                any_cotrans = True
            else:
                s = "N/A"
            line += ' ' + s.rjust(col_width[i])
        if any_cotrans:
            f.write(line+'\n')

        # Third line: main type transporter flux
        key = main_type+"_flux"
        line = key.ljust(maxlen)+" :"
        for i in range(len(resultList)):
            val = resultList[i][key]
            if val == None:
                val = 0.
            line += ' ' + ("% .6g" % val).rjust(col_width[i])
        f.write(line+'\n')

        # Fourth line: main type cotransporter flux
        key = main_type+"_co_flux"
        line = key.ljust(maxlen)+" :"
        any_cotrans = False
        for i in range(len(resultList)):
            val = resultList[i][key]
            if val == None:
                s = 'N/A'
            else:
                s = "% .6g" % val
                any_cotrans = True
            line += ' ' + s.rjust(col_width[i])
        if any_cotrans:
            f.write(line+'\n')

        # Repeat for all other types
        types.remove(main_type)

        for typ in types:
            # Transporter
            line = (typ+"_source").ljust(maxlen)+" :"
            for i in range(len(resultList)):
                line += ' ' + resultList[i][typ].rjust(col_width[i])
            f.write(line+'\n')

            # Cotransporter
            line = (typ+"_cosource").ljust(maxlen)+" :"
            any_cotrans = False
            for i in range(len(resultList)):
                if resultList[i][typ+"_co"] in reactions:
                    s = resultList[i][typ+"_co"]
                    any_cotrans = True
                else:
                    s = "N/A"
                line += ' ' + s.rjust(col_width[i])
            if any_cotrans:
                f.write(line+'\n')

            # Transporter flux
            key = typ+"_flux"
            line = key.ljust(maxlen)+" :"
            for i in range(len(resultList)):
                val = resultList[i][key]
                if val == None:
                    val = 0.
                line += ' ' + ("% .6g" % val).rjust(col_width[i])
            f.write(line+'\n')

            # Cotransporter flux
            key = typ+"_co_flux"
            line = key.ljust(maxlen)+" :"
            any_cotrans = False
            for i in range(len(resultList)):
                val = resultList[i][key]
                if val == None:
                    s = 'N/A'
                else:
                    s = "% .6g" % val
                    any_cotrans = True
                line += ' ' + s.rjust(col_width[i])
            if any_cotrans:
                f.write(line+'\n')

        # Flux distributions ('N/A' in all rows if solution is empty)
        f.write('\n')
        for rea in sorted(reactions):
            index = reactions[rea]
            line = rea.ljust(maxlen)+" :"
            for i in range(len(resultList)):
                try:
                    s = "% .6g" % resultList[i]["solution"][index]
                except IndexError:
                    s = "N/A"
                line += ' ' + s.rjust(col_width[i])
            f.write(line+'\n')


if __name__ == "__main__":
    main()
