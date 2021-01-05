"""
This module defines a number of global constants and helper functions.


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

from builtins import map
from builtins import range
from builtins import object
METANO_VERSION = "1.3.0"

COPYRIGHT_VERSION_STRING = ("metano v" + METANO_VERSION + "\n"
                            "Copyright (C) 2010-2019 Julia Helmecke, Alexander Riemer,\n"
                            "Braunschweig University of Technology,\n"
                            "Dept. of Bioinformatics & Biochemistry.\n\n"
                            "This program comes with ABSOLUTELY NO WARRANTY.\n"
                            "metano is free software distributed under the GNU "
                            "General Public License v3,\nand you are welcome to"
                            " redistribute it under the conditions of that "
                            "license.")

from numpy import isinf

class Verbosity(object):
    """ Levels of verbosity
    """
    NONE           = 0
    ERROR          = 1
    WARNING        = 2
    INFO           = 3
    DEBUG          = 4
    ALL            = 4

class When(object):
    """ Constants for when to perform an action
    """
    NEVER          = 0
    ONCE           = 1
    AUTO           = 2
    ALWAYS         = 3
    ANY            = set((NEVER, ONCE, AUTO, ALWAYS))

def decodeWhen(whenStr, allowedOptions=When.ANY):
    """ decode string specifying when to perform an option
    """
    when = { "never"  : When.NEVER,
             "once"   : When.ONCE,
             "auto"   : When.AUTO,
             "always" : When.ALWAYS }[whenStr.lower()]

    if when not in allowedOptions:
        raise KeyError(whenStr)
    return when


class SolverStatus(object):
    """ Status of optimization
    """
    PREPARED       = -1
    OPTIMAL        =  0
    PRIM_INFEAS    =  1
    DUAL_INFEAS    =  2
    UNBOUNDED      =  3
    TIMEOUT        =  4
    MAXITER        =  5
    NOCONV         =  6
    ERROR          =  7
    UNKNOWN        =  8

def printStatus(status):
    try:
        return { SolverStatus.PREPARED : 'not solved, yet',
                 SolverStatus.OPTIMAL : 'optimal',
                 SolverStatus.PRIM_INFEAS : 'primal infeasible',
                 SolverStatus.DUAL_INFEAS : 'dual infeasible',
                 SolverStatus.UNBOUNDED : 'unbounded',
                 SolverStatus.TIMEOUT : 'time limit exhausted',
                 SolverStatus.MAXITER : 'iteration limit exhausted',
                 SolverStatus.NOCONV : 'not converged',
                 SolverStatus.ERROR : 'error',
                 SolverStatus.UNKNOWN : 'unknown (not converged)' }[status]
    except KeyError:
        return "unexpected value: %r" % status


def cvxpyToSolverStatus(status):
    if isinf(status):
        return SolverStatus.PRIM_INFEAS
    else:
        return SolverStatus.OPTIMAL


def cvxToSolverStatus(status):
    try:
        return { 'optimal' : SolverStatus.OPTIMAL,
                 'primal infeasible' : SolverStatus.PRIM_INFEAS,
                 'dual infeasible' : SolverStatus.DUAL_INFEAS,
                 'unknown' : SolverStatus.UNKNOWN }[status]
    except KeyError:
        return SolverStatus.UNKNOWN


def glpkToSolverStatus(status):
    try:
        return { 'solution is optimal' : SolverStatus.OPTIMAL,
                 'solution is feasible' : SolverStatus.OPTIMAL,
                 'solution is infeasible' : SolverStatus.PRIM_INFEAS,
                 'problem has no feasible solution' : SolverStatus.PRIM_INFEAS,
                 'problem has unbounded solution' : SolverStatus.UNBOUNDED,
                 'solution is undefined' : SolverStatus.UNKNOWN,
                 'e_fault' : SolverStatus.ERROR,
                 'e_objll' : SolverStatus.UNBOUNDED,
                 'e_objul' : SolverStatus.UNBOUNDED,
                 'e_itlim' : SolverStatus.MAXITER,
                 'e_tmlim' : SolverStatus.TIMEOUT,
                 'e_nofeas' : SolverStatus.PRIM_INFEAS,
                 'e_instab' : SolverStatus.ERROR,
                 'e_sing' : SolverStatus.ERROR,
                 'e_noconv' : SolverStatus.NOCONV,
                 'e_nopfs' : SolverStatus.PRIM_INFEAS,
                 'e_nodfs' : SolverStatus.DUAL_INFEAS }[status]
    except KeyError:
        return SolverStatus.UNKNOWN


class ReaFileStruc(object):
    """ Definitions for structure (syntax & some semantics) of reaction file
    """
    delimiter   = ":"    # delimiter separating reaction name and equation
    arrowIrr    = "-->"  # arrow for irreversible reactions
    arrowFlip   = "<--"  # arrow for right-to-left irreversible reactions
    arrowRev    = "<=>"  # arrow for reversible reactions
    commentSign = "#"    # comment sign (text right of sign is ignored)

    class ArrowProperties(object):
        """ Class for storing the properties of a reaction arrow

            Members:
                length     -- length of the arrow string in characters
                reversible -- reversibility flag (False = irreversible reaction)
                flipped    -- flipped flag (True if an irreversible reaction is
                              flipped, i.e. given in reversed order (products on
                              left-hand side, reactants on right-hand side)
        """
        def __init__(self, length=0, reversible=False, flipped=False):
            self.length = length
            self.reversible = reversible
            self.flipped = flipped

        def __str__(self):
            s = "(%u, " % self.length
            if not self.reversible:
                s += "ir"
            s += "reversible, "
            if self.flipped:
                s += "flipped)"
            else:
                s += "normal)"
            return s

    arrow_prop = {arrowIrr : ArrowProperties(len(arrowIrr), False, False),
                  arrowFlip: ArrowProperties(len(arrowFlip), False, True),
                  arrowRev : ArrowProperties(len(arrowRev), True,  False)}


class FbaParam(object):
    """ data structure for passing FBA parameter sets

    Member variables:

    solver           -- name of the solver to be used for the optimization
    maxmin           -- maximize (True) or minimize (False) objective function
    objStr           -- objective function definition (given as string)
    linConstraints   -- list of LinearConstraint objects
    numIter          -- number of iterations in case of nonlinear optimization
    nlc              -- list of nonlinear constraints (lambda functions)
    nlc_grad         -- gradients of nonlinear constraints
    """
    DEFAULT_NUMITER = 10

    def hasEqualityConstraints(self):
        for c in self.linConstraints:
            if c.isEq:
                return True
        return False

    def __init__(self, solver="default", maxmin=True, objStr="", numIter=-1):
        self.solver  = solver
        self.maxmin  = maxmin
        self.objStr  = objStr
        self.numIter = -1
        self.linConstraints = []
        self.nlc = []
        self.nlc_grad = []

    def setLinConstraints(self, lc):
        self.linConstraints = lc

    def setNLConstraints(self, nlc, nlc_grad):
        self.nlc = nlc
        self.nlc_grad = nlc_grad

    def __repr__(self):
        return "FbaParam(%r, %r, %r, %r)" % (self.solver, self.maxmin,
                                             self.objStr, self.numIter)


def typename(x):
    """ return the most informative type name available for x
    """
    typeOfX = type(x).__name__
    if typeOfX == "instance":
        return x.__class__.__name__
    return typeOfX


def printHistogram(histogram, minval, maxval, revert=False):
    numBins = len(histogram)
    if not numBins:
        return
    if revert:
        tmp = -minval
        minval = -maxval
        maxval = tmp

    print ("\nmin value (0): %.12g, max value (%u): %.12g\nHistogram:\n" %
           (minval, numBins-1, maxval))

    binWidth = (maxval-minval)/float(numBins)
    if revert:
        # Print histogram negated and in reverse order
        for i in range(numBins):
            print(("%.12g" % (minval + float(i)*binWidth)).ljust(16), end=' ')
            print(-histogram[-i-1])
    else:
        for i in range(numBins):
            print(("%.12g" % (minval + float(i)*binWidth)).ljust(16), end=' ')
            print(histogram[i])
    print(("%.12g" % maxval).ljust(16))


def exportMatrixAsDotM(matrix, filename):
    """ export a matrix as .m file (for Matlab/Octave)
    """
    with open(filename, 'w') as f:
        exportMatrixAsDotMByHandle(matrix, f)


def exportMatrixAsDotMByHandle(matrix, f):
    """ export a matrix as .m file (for Matlab/Octave) to the given file handle
        (must be open for writing)
    """
    nRows = len(matrix)
    if not nRows:
        return      # nothing to do

    f.write("S = [")
    i = 0
    for row in matrix:
        f.write(' '.join(map(repr, row)))
        i += 1
        if i >= nRows:
            f.write(']')
        f.write(";\n")


def makeCoefNameTerm(s, coef, name):
    """ generate a string representation of the (coef, name) pair with an
        appropriate sign connecting the term to preceding terms

    Keyword arguments:

    s    -- string containing the preceding terms (may be empty)
    coef -- coefficient (positive or negative floating-point value)
    name -- name string to be displayed next to coefficient
    """
    # Suppress terms with a coefficient of zero, and omit coefficient 1.
    if coef == 0.:
        return ""
    absCoef = abs(coef)
    if absCoef == 1.:
        absCoefTerm = "%s" % name
    else:
        absCoefTerm = "%r %s" % (absCoef, name)

    # First coefficient is marked with minus sign or nothing, others with " + "
    # or " - "
    if s == "":
        if coef < 0:
            sign = "-"
        else:
            sign = ""
    else:
        if coef < 0:
            sign = " - "
        else:
            sign = " + "

    return sign + absCoefTerm


def displayPermutedMatrix(matrix, permutation, names):
    """ return a list of lines for displaying the given sparsely populated
        matrix

    Output is in form of linear combinations with the coefficients coming from
    the matrix and the variables coming from the list of names. Example:

      displayPermutedMatrix([[0, 1, -1]], [1, 0, 2], ["A", "B", "C"]) yields
      ["1 A - 1 C"]

    Keyword arguments:

    matrix      -- matrix
    permutation -- permutation[i] is the entry in names for matrix column i
    names       -- list of reaction names
    """
    nCols = len(permutation)
    lines = []
    for row in matrix:
        s = ""
        for i in range(nCols):
            coef = row[i]
            if coef:
                s += makeCoefNameTerm(s, coef, names[permutation[i]])
        lines.append(s)
    return lines


def makeUnique(strings):
    """ make a list of strings unique by appending arbitrary numbers,
        e.g. if string "abc" appears twice in the list, the occurrences are
        replaced with "abc_001" and "abc_002"
    """
    result = list(strings)  # copy list
    nStrings = len(strings)
    sDict = {}  # dict { string : list of indices }
    for i in range(nStrings):
        name = strings[i]
        if name in sDict:
            sDict[name].append(i)
        else:
            sDict[name] = [i]

    if len(sDict) != nStrings:
        for name in sDict:
            indexList = sDict[name]
            if len(indexList) > 1:
                for i in range(len(indexList)):
                    # simply add number to name
                    result[indexList[i]] += "_%03u" % (i+1)
    return result


def addToDictPlus(d, key, value):
    """ if key doesn't exist, add (key : value) pair to d, otherwise add value
        to d[key]
    """
    if key in d:
        d[key] += value
    else:
        d[key] = value


def roundAndCut(number, ndigits=8, threshold=1E-8):
    """ round a number to ndigits digits and set it to zero if its absolute
        value is below threshold
    """
    if abs(number) < threshold:
        return 0.
    else:
        return round(number, ndigits)

def isZeroVector(vec):
    for elem in vec:
        if elem != 0.:
            return False
    return True

def padNumber(numStr):
    """ pad a number given as a string representation by adding a space to the
        left if the number is not negative - string must not be empty
    """
    if numStr[0] == '-':
        return numStr
    else:
        return ' ' + numStr
