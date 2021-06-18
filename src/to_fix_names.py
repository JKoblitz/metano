#!/usr/bin/env python

""" Script for fixing a reaction file

This script constructs a ReactionParser and calls its fixfile() method,
which replaces spaces in reaction or metabolite names with underscores.


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
from metano.defines import Verbosity, COPYRIGHT_VERSION_STRING
import os


def main():
    # Parse command line

    usage = "Usage: %prog [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-i", "--input", dest="inputFile", help="use the "
                      "given reaction FILE as input", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile", help="write output "
                      "(reaction file) to FILE (must not be the same as input "
                      "file)", metavar="FILE")

    options, _ = parser.parse_args()
    parser.check_required("-i")
    parser.check_required("-o")

    if os.path.exists(options.outputFile):
        if os.path.samefile(options.outputFile, options.inputFile):
            print("Error: Output file must not be the same as input file!")
            exit()

    # Call reaction parser for fixing the reaction file

    rparser = ReactionParser()
    try:
        rparser.fixfile(options.inputFile, options.outputFile)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s or write file "
               "%s:" % (os.path.normpath(options.inputFile),
                        os.path.normpath(options.outputFile)))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in reaction file %s:" %
               os.path.normpath(options.inputFile))
        print(strerror)
        exit()

    # Print warning and info messages
    print('\n'.join(rparser.getMessages(Verbosity.ALL)))


if __name__ == "__main__":
    main()
