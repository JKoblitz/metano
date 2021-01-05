"""
This module defines an iterator that returns a given number of random vectors.


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

from builtins import object
from numpy.random import uniform


class StartPointIterator(object):
    """ Iterator for generating random start points, i.e. vectors with dimension
        n and entries from a uniform distribution in [lb, ub)
    """
    def __init__(self, dimension, numIterations):
        self.dim = dimension
        self.numIter = numIterations
        self.lb = -1.
        self.ub = 1.

    def setDimension(self, dimension):
        self.dim = dimension

    def setRange(self, lb=None, ub=None):
        if lb != None:
            self.lb = lb
        if ub != None:
            self.ub = ub

    def reset(self, numIterations):
        self.numIter = numIterations

    def __iter__(self):
        return self

    def __next__(self):
        if self.numIter == 0:
            raise StopIteration
        self.numIter -= 1
        return uniform(self.lb, self.ub, self.dim)
