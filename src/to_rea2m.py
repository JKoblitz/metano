#!/usr/bin/env python

""" Script for converting a reaction file to a stoichiometric matrix

This script constructs a ReactionParser and calls its exportMatrixAsDotM()
method, which exports the stoichiometric matrix as a .m file.


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

from fba import OptionParser
from metano.metabolicmodel import MetabolicModel
from metano.defines import exportMatrixAsDotM, COPYRIGHT_VERSION_STRING
import os


def main():
    # Parse command line

    usage = "Usage: %prog [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile", help="use the "
                      "given reaction FILE as input", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile", help="write the "
                      "stoichiometric matrix to FILE (should have .m extension)",
                      metavar="FILE")

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-o")

    # Parse reaction file

    model = MetabolicModel()
    try:
        model.addReactionsFromFile(options.reactionFile)
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

    # Export matrix

    try:
        exportMatrixAsDotM(model.getStoichiometricMatrix(), options.outputFile)
    except IOError as strerror:
        print ("An error occurred while trying to write file %s:" %
               os.path.basename(options.outputFile))
        print(strerror)
        exit()


if __name__ == "__main__":
    main()
