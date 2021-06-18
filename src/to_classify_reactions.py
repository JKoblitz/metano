#!/usr/bin/env python
"""

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
from metano.metabolicmodel import MetabolicModel
from metano.defines import COPYRIGHT_VERSION_STRING
import os
import sys


def main():

    # 1. Parse command line

    usage = "Usage: %prog <filename> [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-c", "--compartment-suffixes",
                      dest="compartmentSuffixes", help="STRING of compartment "
                      "suffixes (separated by whitespace)", metavar="STRING")

    options, args = parser.parse_args()

    if len(args) < 1:
        print("Error: No filename given.")
        print ("Usage is\n    " + os.path.basename(sys.argv[0]) + " <filename> "
               "[options]")
        exit()

    filename = args[0]

    # 2. Parse reaction file

    model = MetabolicModel()
    try:
        model.addReactionsFromFile(filename)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(filename))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in reaction file %s:" %
               os.path.basename(filename))
        print(strerror)
        exit()

    # 3. Identify exchange reactions

    exchangeReactions = [reaction for reaction in model if len(reaction) == 1]
    nExchange = len(exchangeReactions)
    print("Exchange reactions:", nExchange)
    print("\n".join([r.name for r in exchangeReactions
                     if "exchange" not in r.name]))

    # 4. Identify transport reactions

    if not options.compartmentSuffixes:
        print ("Warning: No compartment suffixes given. Unable to identify "
               "transport reactions.")


if __name__ == "__main__":
    main()
