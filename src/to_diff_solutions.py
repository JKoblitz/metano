#!/usr/bin/env python

""" This script computes the absolute difference between two flux distributions
    given as FBA solution files.

Usage:

    to_diff_solutions.py <file1> <file2> [-o output-file] [-c cutoff]

The output is written to the console if no output file is given.
Output format is as follows (without header):
<NAME>    <DIFFERENCE>

Reactions are sorted in order of decreasing flux difference. Reactions in which
the fluxes differ by less than cutoff are omitted.


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

from builtins import map
from builtins import range
from past.utils import old_div
from metano.fba import OptionParser
from metano.metabolicflux import MetabolicFlux
from metano.defines import COPYRIGHT_VERSION_STRING
import os
import sys


def main():
    # 1. Parse command line

    usage = "Usage: %prog <file1> <file2> [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-o", "--output", dest="outputFile", help="output FILE "
                      "(output is written to console if this is missing)",
                      metavar="FILE")
    parser.add_option("-c", "--cutoff", dest="cutoff", help="Cutoff below which"
                      " difference is considered zero; default is 1E-10")
    parser.set_defaults(cutoff="1E-10")

    options, args = parser.parse_args()

    if len(args) < 2:
        print("Error: Need two solution files.")
        print ("Usage is\n    " + os.path.basename(sys.argv[0]) + " <file1> "
               "<file2> [options]")
        exit()

    try:
        cutoff = float(options.cutoff)
    except ValueError:
        print("Error: Invalid floating point value for cutoff.")
        exit()
    if cutoff < 0.:
        print ("Warning: Cutoff is less than zero. Setting cutoff to zero, "
               "i.e. no cutoff.")
        cutoff = 0.

    basename = list(map(os.path.basename, args))
    # Use full path for identifying files if basenames are identical
    if basename[0] == basename[1]:
        basename = args

    # 2. Parse solution files and compute absolute differences

    solution = MetabolicFlux(), MetabolicFlux()
    for i in range(2):
        # Get pair of MetabolicFlux objects
        try:
            solution[i].readFromFile(args[i])
        except IOError as strerror:
            print ("An error occurred while trying to read file %s:" %
                   os.path.basename(basename[i]))
            print(strerror)
            exit()
        except SyntaxError as strerror:
            print("An error occurred parsing file %s:" % basename[i])
            print(strerror)
            exit()

    diffSolution = solution[0].absDiff(solution[1])

    # 3. Generate output

    nReactions = len(diffSolution)
    if nReactions == 0:
        exit()  # Nothing to do
    if len(diffSolution) == len(solution[0]) == len(solution[1]):
        print("Both solutions have the same %u reactions." % nReactions)
    else:
        n = list(map(len, solution))
        print ("The two solutions have %u reactions in common.\n(%s: %u, "
               "%s: %u)" % (n[0]+n[1]-nReactions, basename[0], n[0],
                            basename[1], n[1]))

    sumDiff = sum(diffSolution.fluxDict.values())
    sumSqDiff = sum(x*x for x in list(diffSolution.fluxDict.values()))
    print("Summed absolute difference: %s (mean: %s)" % (sumDiff,
                                                         old_div(sumDiff, nReactions)))
    print("Summed squared difference: %s (mean: %s)" % (sumSqDiff,
                                                        old_div(sumSqDiff, nReactions)))

    try:
        if options.outputFile is None:
            print()
            written = write_output(sys.stdout, diffSolution, solution, basename,
                                   cutoff, False)
        else:
            with open(options.outputFile, 'w') as f:
                written = write_output(f, diffSolution, solution, basename,
                                       cutoff, True)
    except IOError as strerror:
        print("An error occurred while trying to write output:")
        print(strerror)
        exit()

    if not written:
        print ("All differences are below the cutoff (%g). Nothing written." %
               cutoff)


def write_output(file_obj, diffSolution, solution, filenames, cutoff, showOrig):
    """ Write output to the given file object (expected to be open for writing)

    This function writes a list of "<reaction> <difference>" pairs to the given
    file object, sorted in order of decreasing difference. Differences below
    cutoff are not reported. If cutoff is zero, all differences are reported.
    If showOrig is True, the original flux values are listed in addition to the
    absolute differences.

    Keyword arguments:
    file_obj     -- the file object (open for writing; may also be sys.stdout)
    diffSolution -- MetabolicFlux of absolute flux differences
    solution     -- pair (or list) of solutions (flux distributions),
                    as MetabolicFlux
    filenames    -- pair (or list) of solution filenames, in original order
    cutoff       -- cutoff below which differences are not reported
    showOrig     -- if True, also show original values, else only the difference

    Returns:
    True if any difference is above cutoff, False else
    """

    # Get length of longest reaction name and length of filenames for formatting
    maxlen = len(max(list(diffSolution.fluxDict.keys()), key=len))
    len0 = max(len(filenames[0]), 12)
    len1 = max(len(filenames[1]), 12)

    first = True
    # Sort dictionary by value in decreasing order
    for rea in sorted(diffSolution.fluxDict, key=diffSolution.fluxDict.get,
                      reverse=True):
        value = diffSolution[rea]
        if value < cutoff:
            break

        if first:
            s = "REACTION".ljust(maxlen)+"DIFFERENCE".rjust(12)
            if showOrig:
                s += ('  '+filenames[0].rjust(len0)+'  ' +
                      filenames[1].rjust(len1))
            file_obj.write(s+'\n')
            first = False

        s = rea.ljust(maxlen) + ("%.6g" % value).rjust(12)
        in2 = rea in solution[1]
        if rea in solution[0]:
            if showOrig:
                s += '  '+("% .6g" % solution[0][rea]).rjust(len0)
            notPresent = -1
        else:
            notPresent = 0
            if showOrig and in2:
                s += ' '*(len0+2)
        if in2:
            if showOrig:
                s += '  '+("% .6g" % solution[1][rea]).rjust(len1)
        else:
            notPresent = 1

        if notPresent >= 0:
            s += " # not present in " + filenames[notPresent]
        file_obj.write(s+'\n')

    return not first


if __name__ == "__main__":
    main()
