#!/usr/bin/env python

"""
Script for adding bounds imposed by the irreversibility of some reactions to
a metano scenario file

This script reads a scenario file and a reaction file and writes a scenario
file. While the two scenario files can be the same, it is advisable to keep
the original scenario file, which is otherwise overwritten.

The central function is performed by
MetabolicModel.writeBoundsToParamFileHandle().


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
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.defines import COPYRIGHT_VERSION_STRING
import os


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-p", "--parameters", dest="inputFile",
                      help="use the given scenario FILE", metavar="FILE")
    parser.add_option("-r", "--reactions", dest="reactionFile", help="use "
                      "reversibilities from reaction FILE", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile", help="write output "
                      "(scenario file) to FILE", metavar="FILE")

    options, _ = parser.parse_args()
    parser.check_required("-p")
    parser.check_required("-r")
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
        maxmin, obj_name, solver, numIter, lb, ub = \
            pparser.parse(options.inputFile)
        # Set flux bounds in model (for parameter check)
        model.setFiniteBounds(lb, ub, True, model_messages)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.inputFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in scenario file %s:" %
               os.path.basename(options.inputFile))
        print(strerror)
        exit()
    except ValueError as strerror:
        print(strerror)
        exit()

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print('\n'.join(msgs))

    # 4. Write new scenario file

    try:
        with open(options.outputFile, 'w') as f:
            str_maxmin = {True: "max", False: "min"}[maxmin]
            f.write("OBJ %s %s\n" % (str_maxmin, obj_name))
            if solver != "":
                f.write("SOLVER %s\n" % solver)
            if numIter >= 0:
                f.write("NUM_ITER %u\n" % numIter)
            f.write("\n")
            model.writeBoundsToParamFileHandle(f)
            f.write("\n")
    except IOError as strerror:
        print ("An error occurred while trying to write file %s:" %
               os.path.basename(options.outputFile))
        print(strerror)
        exit()


if __name__ == "__main__":
    main()
