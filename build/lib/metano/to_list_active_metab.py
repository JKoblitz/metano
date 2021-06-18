#!/usr/bin/env python

""" Script for determining the metabolites in active reactions

This script parses an FBA solution file and the underlying reaction file and
prints all metabolites present in reactions with a non-zero flux along with
the number of such reactions in which they appear. It also lists the active
reactions.


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

from metano.fba import OptionParser
from metano.reactionparser import ReactionParser
from metano.metabolicmodel import MetabolicModel
from metano.metabolicflux import MetabolicFlux
from metano.defines import COPYRIGHT_VERSION_STRING
import os


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="use the given reaction FILE", metavar="FILE")
    parser.add_option("-s", "--solution", dest="solutionFile",
                      help="parse the given FBA solution FILE", metavar="FILE")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float",
                      help="threshold below which a flux is to be considered "
                      "zero; default: 1E-10)", metavar="VALUE")
    parser.set_defaults(cutoff=1E-10)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-s")

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

    # 3. Parse solution file

    solution = MetabolicFlux()
    try:
        solution.readFromFile(options.solutionFile)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.solutionFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in solution file %s:" %
               os.path.basename(options.solutionFile))
        print(strerror)
        exit()

    # 4. Count metabolites' participation in active reactions

    model_messages = []
    activeReactions, metabolitesSolution = \
        model.getActiveReactions(solution.getVecOrderedByModel(model),
                                 options.cutoff, model_messages)

    # Show warning and info messages (if any)
    msgs = rparser.getMessages() + [x[1] for x in model_messages]
    if msgs:
        print('\n'.join(msgs)+'\n')

    # 5. Output results

    print("Active metabolites:\n")
    for met in sorted(metabolitesSolution, key=metabolitesSolution.get,
                      reverse=True):
        print(met, metabolitesSolution[met])

    print("\nActive reactions:\n")
    for reactionName in activeReactions:
        print(reactionName)

    print ("\n%u metabolites in %u active reactions." %
           (len(metabolitesSolution), len(activeReactions)))


if __name__ == "__main__":
    main()
