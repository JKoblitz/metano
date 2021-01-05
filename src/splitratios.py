#!/usr/bin/env python

""" This script computes the split ratios for the fluxes entering and leaving
    each metabolite (as percentages).

Usage:

    splitratios.py <solution-file> -r reaction-file [-o output-file]
                   [-m metabolite] [-t threshold] [-u] [-c cutoff] [-f] [-a]


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
import os
import sys
import csv
from metano.fba import OptionParser
from metano.reactionparser import ReactionParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.defines import COPYRIGHT_VERSION_STRING


def main():
    # 1. Parse command line

    usage = "Usage: %prog <solution-file> [options]"
    version = "Split-ratio analysis\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile", help="use "
                      "metabolic model given by reaction FILE", metavar="FILE")
    parser.add_option("-m", "--metabolite", dest="metabolite", help="NAME of "
                      "the single metabolite to be analyzed (optional)",
                      metavar="NAME")
    parser.add_option("-o", "--output", dest="outputFile", help="perform "
                      "analysis for all metabolites and write output to FILE",
                      metavar="FILE")
    parser.add_option("-t", "--omit-fluxes-below", dest="threshold",
                      type="float", help="omit metabolites in output whose flux"
                      " is below the given THRESHOLD (only with -o)",
                      metavar="THRESHOLD")
    parser.add_option("-u", "--omit-unbranched", dest="omitUnbranched",
                      action="store_true", help="if set, metabolites with "
                      "unbranched flux are omitted (only with -o)")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float", help="use "
                      "the given cutoff VALUE (default 0)", metavar="VALUE")
    parser.add_option("-f", "--cutoff-is-absolute", dest="cutoffIsAbsolute",
                      action="store_true", help="if set, cutoff is in absolute "
                      "flux units rather than a value between 0 and 1")
    parser.add_option("-a", "--list-all", dest="listAll", action="store_true",
                      help="if set, cutoff is ignored, and even zero fluxes are"
                      " shown (per definition as outgoing)")
    parser.set_defaults(threshold=-1., omitUnBranched=False, cutoff=0.,
                        cutoffIsAbsolute=False, listAll=False)

    options, args = parser.parse_args()
    parser.check_required('-r')

    if not options.metabolite and not options.outputFile:
        print("Neither metabolite nor output file given. Nothing to do.")
        exit()

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

    # 3. Parse reaction file

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

    # 4. Compute split ratios

    if options.outputFile:
        # Compute split ratios for all metabolites
        splitRatios = solution.computeAllSplitRatios(model, options.cutoff,
                                                     options.cutoffIsAbsolute, options.listAll)
        try:
            writeSplitRatiosToFile(options.outputFile, splitRatios,
                                   options.threshold, options.omitUnbranched)
        except IOError as strerror:
            print ("Unable to write to file %s:" %
                   os.path.basename(options.outputFile))
            print(strerror)

        # Print statistics
        nMetabolites = len(splitRatios)
        print ("The model contains %u reactions and %u metabolites." %
               (len(model), nMetabolites))
        EPSILON = 1e-6
        nZero = nOne = nBelowEps = nDiff = 0
        maxDiff = 0.
        maxdMet = ""

        for met in splitRatios:
            outRatios, inRatios = splitRatios[met]
            outSum = sum(outRatios[rea][1] for rea in outRatios)
            diffInOut = (abs(outSum - sum(inRatios[rea][1]
                                          for rea in inRatios)))
            if diffInOut > EPSILON:
                nDiff += 1
                if diffInOut > maxDiff:
                    maxDiff = diffInOut
                    maxdMet = met
            if len(outRatios) == 0 or outSum == 0.:
                nZero += 1
            else:
                if len(outRatios) == 1 and len(inRatios) == 1:
                    nOne += 1
                elif options.listAll:
                    # Actually count all outgoing and incoming fluxes > 0
                    nOut = sum(
                        1 for rea in outRatios if outRatios[rea][1] > 0.)
                    nIn = sum(1 for rea in inRatios if inRatios[rea][1] > 0.)
                    if nOut == nIn == 1:
                        nOne += 1
                if sum(outRatios[rea][1] for rea in outRatios) < EPSILON:
                    nBelowEps += 1

        print ("%u metabolites have a diff. betw. incoming and outgoing "
               "fluxes above %g.\nThe highest diff. betw. incoming and "
               "outgoing metabolite flux is %g (%s)." % (nDiff, EPSILON,
                                                         maxDiff, maxdMet))
        if options.cutoff == 0.:
            print ("%u metabolites (%.2g%%) have flux zero." %
                   (nZero, (float(nZero)/float(nMetabolites))*100.))
        else:
            if options.cutoffIsAbsolute:
                sCutoff = "%g" % options.cutoff
            else:
                sCutoff = "%g%%" % (options.cutoff * 100.)
            print ("%u metabolites (%.2g%%) have no single outgoing flux above "
                   "%s." % (nZero, (float(nZero)/float(nMetabolites))*100.,
                            sCutoff))
        if nBelowEps:
            print ("Another %u metabolites (%.2g%%) have a flux below %g." %
                   (nBelowEps, ((float(nBelowEps)/float(nMetabolites))*100.),
                    EPSILON))
        print ("%u metabolites (%.2g%%) lie on an unbranched path." %
               (nOne, ((float(nOne)/float(nMetabolites))*100.)))
        print()

        # Extract data for queried metabolite
        if options.metabolite:
            try:
                outRatios, inRatios = splitRatios[options.metabolite]
            except KeyError:
                print ("Error: Metabolite %r does not exist." %
                       options.metabolite)
                exit()

    elif options.metabolite:
        # Only compute split ratios for the given metabolite
        outRatios, inRatios = solution.computeSplitRatios(options.metabolite,
                                                          model, options.cutoff, options.cutoffIsAbsolute, options.listAll)

    if options.metabolite:
        print("Split ratios for %r" % options.metabolite)
        maxlen = max(len(max(outRatios, key=len)), len(max(inRatios, key=len)))
        fluxSum = sum(outRatios[rea][1] for rea in outRatios)
        print("- outgoing (sum %g) -" % fluxSum)
        for rea in sorted(outRatios, key=lambda x: outRatios[x][0],
                          reverse=True):
            ratio, flux = outRatios[rea]
            print (rea.ljust(maxlen) + " " + ("%.4g%%" % (ratio*100.)).rjust(10)
                   + " " + ("%.6g" % flux).rjust(12))
        fluxSum = sum(inRatios[rea][1] for rea in inRatios)
        print("- incoming (sum %g) -" % fluxSum)
        for rea in sorted(inRatios, key=lambda x: inRatios[x][0],
                          reverse=True):
            ratio, flux = inRatios[rea]
            print (rea.ljust(maxlen) + " " + ("%.4g%%" % (ratio*100.)).rjust(10)
                   + " " + ("%.6g" % flux).rjust(12))


def writeSplitRatiosToFile(filename, splitRatios, threshold, omitUnbranched):
    with open(filename, 'w') as f:
        if filename.split(".")[-1] == "csv":
            writeSplitRatiosToCSVHandle(f, splitRatios, threshold,
                                        omitUnbranched)
        else:
            writeSplitRatiosToFileHandle(f, splitRatios, threshold,
                                         omitUnbranched)


def writeSplitRatiosToCSVHandle(f, splitRatios, threshold, omitUnbranched):
    # sorted[Out/In]Keys stores the sorted keys for the outRatios or inRatios
    #  dict of each metabolite (in descending order of flux)
    sortedOutKeys, sortedInKeys = {}, {}
    # maxwidth stores maximum width of each column
    metaboliteList = sorted(splitRatios)
    csvwriter = csv.writer(f, delimiter=";", quotechar='"',
                           quoting=csv.QUOTE_MINIMAL)
    # Write Header to CSV file
    csvwriter.writerow(["NAME", "DIRECTION", "REACTION", "PERCENT", "FLUX"])

    for met in sorted(metaboliteList):
        outRatios, inRatios = splitRatios[met]

        if omitUnbranched and len(outRatios) <= 1 and len(inRatios) <= 1:
            # Omit metabolite in output if it has an unbranched flux
            continue

        if not outRatios:
            csvwriter.writerow([met, "outgoing", "", 0, 0])
        if not inRatios:
            csvwriter.writerow([met, "incoming", "", 0, 0])
        for o in sorted(outRatios, key=outRatios.get, reverse=True):
            csvwriter.writerow(
                [met, "outgoing", o, outRatios[o][0], outRatios[o][1]])
        for i in sorted(inRatios, key=inRatios.get, reverse=True):
            csvwriter.writerow(
                [met, "incoming", i, inRatios[i][0], inRatios[i][1]])


def writeSplitRatiosToFileHandle(f, splitRatios, threshold, omitUnbranched):
    # sorted[Out/In]Keys stores the sorted keys for the outRatios or inRatios
    #  dict of each metabolite (in descending order of flux)
    sortedOutKeys, sortedInKeys = {}, {}
    # maxwidth stores maximum width of each column
    maxwidth = {}
    metaboliteList = sorted(splitRatios)

    for met in metaboliteList:
        outRatios, inRatios = splitRatios[met]
        sortedOutKeys[met] = sorted(outRatios, key=lambda x: outRatios[x][0],
                                    reverse=True)
        sortedInKeys[met] = sorted(inRatios, key=lambda x: inRatios[x][0],
                                   reverse=True)
        for i in range(len(sortedOutKeys[met])):
            try:
                maxwidth[i] = max(maxwidth[i], len(sortedOutKeys[met][i]))
            except KeyError:
                maxwidth[i] = max(len(sortedOutKeys[met][i]), 12)
        for i in range(len(sortedInKeys[met])):
            try:
                maxwidth[i] = max(maxwidth[i], len(sortedInKeys[met][i]))
            except KeyError:
                maxwidth[i] = max(len(sortedInKeys[met][i]), 12)

    for met in metaboliteList:
        outRatios, inRatios = splitRatios[met]

        if omitUnbranched and len(outRatios) <= 1 and len(inRatios) <= 1:
            # Omit metabolite in output if it has an unbranched flux
            continue

        # 1st line: <metabolite name>
        # 2nd line: -outgoing <flux>
        # 3rd line: <rea1> <rea2> ... <reaN>
        # 4th line:  <%1>   <%2>  ...  <%N>
        # 5th line:  <v1>   <v2>  ...  <vN>
        # 6th line: -outgoing <flux>
        # 7th line: <rea1> <rea2> ... <reaN>
        # 8th line:  <%1>   <%2>  ...  <%N>
        # 9th line:  <v1>   <v2>  ...  <vN>
        lines = ["", "", ""]
        fluxSum = 0.
        for i in range(len(sortedOutKeys[met])):
            rea = sortedOutKeys[met][i]
            ratio, flux = outRatios[rea]
            fluxSum += flux
            lines[0] += " " + rea.rjust(maxwidth[i])
            lines[1] += " " + ("%.4g%%" % (ratio*100.)).rjust(maxwidth[i])
            lines[2] += " " + ("%.6g" % flux).rjust(maxwidth[i])

        if fluxSum < threshold:
            # Omit metabolite in output if its total flux is below threshold
            continue

        s = "%s\n-outgoing %g\n" % (met, fluxSum)
        for line in lines:
            if line:
                s += line + '\n'

        lines = ["", "", ""]
        fluxSum = 0.
        for i in range(len(sortedInKeys[met])):
            rea = sortedInKeys[met][i]
            ratio, flux = inRatios[rea]
            fluxSum += flux
            lines[0] += " " + rea.rjust(maxwidth[i])
            lines[1] += " " + ("%.4g%%" % (ratio*100.)).rjust(maxwidth[i])
            lines[2] += " " + ("%.6g" % flux).rjust(maxwidth[i])

        s += "-incoming %g\n" % fluxSum
        for line in lines:
            if line:
                s += line + '\n'
        f.write(s)


if __name__ == "__main__":
    main()
