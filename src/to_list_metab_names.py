#!/usr/bin/env python
""" Script for listing all metabolite names

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

from metano.metabolicmodel import MetabolicModel
import sys
import os

try:
    filename = sys.argv[1]
except IndexError:
    print("Error: No filename given.")
    print("Usage is\n    " + os.path.basename(sys.argv[0]), "<reaction file>")
    exit()

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
print("\n".join(sorted(model.getMetaboliteNames())))
