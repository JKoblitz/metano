"""
This file defines a set of functions and classes for cycle analysis
according to
[1] RR Vallabhajosyula, V Chickarmane, HM Sauro: Conservation analysis of
    large biochemical networks. Bioinformatics 22(3): 346-353 (2006)

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
from __future__ import absolute_import

from builtins import range
from builtins import object
from metano.defines import roundAndCut, makeCoefNameTerm, displayPermutedMatrix
from numpy import array, zeros, empty, eye, nonzero, vectorize
try:
    import pygsl
    import pygsl.permutation
    import pygsl.linalg
    _pygsl_avail = True
except ImportError:
    _pygsl_avail = False


def _myround(x): return roundAndCut(x, 8)


_vmyround = vectorize(_myround)


class FluxCycle(object):
    """ class for representing flux cycles

    Member variables:

    lhs     - matrix
    rhs     - vector or list of strings
    indices - array of indices linking the columns of lhs to the full list of
              reactions (indices[i] is the index of column i in the original
              list

    A flux cycle is represented by a linear system of equations, which is stored
    in the variables lhs and rhs for the left-hand and right-hand sides.
    """

    def __init__(self, lhs=None, rhs=None, indices=[]):
        """ construct a new FluxCycle object
        """
        self.lhs = lhs
        self.rhs = rhs
        self.indices = indices

    def display(self, reactionNames):
        """ return a list of lines for displaying the linear system of equations
        """
        if self.lhs is None or self.rhs is None:
            return []

        # First convert left-hand side to (coefficient, name) pairs joined by
        # plus or minus signs
        lines = displayPermutedMatrix(self.lhs, self.indices, reactionNames)

        try:
            float(self.rhs[0])
            isString = False
        except ValueError:
            isString = True
        except IndexError:
            return lines

        # Pad lines to same length and add right-hand side
        maxlen = len(max(lines, key=len))
        nRows = len(self.lhs)
        if isString:
            for i in range(nRows):
                if self.rhs[i] == "":
                    lines[i] = lines[i].ljust(maxlen) + " = 0"
                else:
                    lines[i] = "%s = %s" % (
                        lines[i].ljust(maxlen), self.rhs[i])
        else:
            for i in range(nRows):
                lines[i] = "%s = %g" % (lines[i].ljust(maxlen), self.rhs[i])

        return lines


def gaussJordan(M, epsilon=1E-8):
    """ transform the given matrix to reduced row-echelon form using Gauss-
        Jordan elimination and read off the rank of M. The operation is
        performed in place, i.e. M is overwritten with the result.
    """
    lead = 0
    nRows, nCols = M.shape
    rank = nRows  # assume full rank

    for i in range(nRows):
        # Find leading coefficient in row i
        while lead < nCols and abs(M[i, lead]) < epsilon:
            # Look for non-zero pivot element in other rows
            j = i+1
            while j < nRows and abs(M[j, lead]) < epsilon:
                j += 1

            if j < nRows:
                # Swap row j with row i
                temp = M[i, :].copy()
                M[i, :] = M[j, :]
                M[j, :] = temp
                break

            #  If all elements in column are zero, proceed to next column
            lead += 1

        if lead >= nCols:
            rank = i  # rank is index of first all-zero row
            break

        M[i, :] /= M[i, lead]  # Divide row i by leading coefficient

        # Subtract scaled row i from all rows below and above to set all
        # M[j, lead] to 0
        for j in range(i):
            M[j, :] -= M[j, lead] * M[i, :]
        for j in range(i+1, nRows):
            M[j, :] -= M[j, lead] * M[i, :]

        lead += 1

    return rank


def ref2rref(M, epsilon=1E-8):
    """ convert a matrix M that is already in row-echelon form to reduced
        row-echelon form and read off the rank of M. The operation is performed
        in place, i.e. M is overwritten with the result.
    """
    lead = 0
    nRows, nCols = M.shape
    rank = nRows  # assume full rank

    for i in range(nRows):
        # Find leading coefficient in row i
        while lead < nCols and abs(M[i, lead]) < epsilon:
            lead += 1

        if lead >= nCols:
            rank = i  # rank is index of first all-zero row
            break

        M[i, :] /= M[i, lead]  # Divide row i by leading coefficient

        # Subtract scaled row i from all rows above to set all M[j, lead] to 0
        for j in range(i):
            M[j, :] -= M[j, lead] * M[i, :]

        lead += 1

    return rank


def qrpt_decomp2(A):
    """ compute the QR factorization with column pivoting of the matrix A
        (wrapper for the corresponding GSL function)

    Keyword arguments:

    A     -- input matrix

    Returns tuple (p, q, r) with

    p     -- permutation P as list (column index in R -> original index in A)
    q     -- orthogonal matrix Q
    r     -- upper triangular matrix R with lower (all-zero) part removed

    with  QRP' = A
    """
    m, n = A.shape
    dtype = A.dtype
    q = zeros((m, m), dtype)
    r = zeros((m, n), dtype)
    tau = pygsl.permutation.Permutation(min(m, n))
    p = pygsl.permutation.Permutation(n)
    sig = pygsl.permutation.Permutation(n)

    # Perform QR decomposition with column pivoting on matrix A
    pygsl.linalg._gslwrap.gsl_linalg_QRPT_decomp2(A, q, r, tau, p, sig)

    # Transform R to reduced row-echelon form and read off rank
    rank = ref2rref(r)

    return p.toarray(), q, r[:rank, :]


class CycleAnalyzer(object):
    """ class for cycle analysis based on [1]

    Member variables:

    gamma       -- matrix of linear dependencies as defined in [1]
    gIndices    -- indices of columns of gamma in original stoichiometric matrix
    r           -- reduced row-echelon form of permuted stoichiometric matrix
    rIndices    -- indices of columns of r in original stoichiometric matrix
    nCycleReacs -- number of reactions that are involved in at least one cycle
                   rIndices[:nCycleReacs] are the indices of these
    """

    def __init__(self):
        """ construct a new CycleAnalyzer object
        """
        self.gamma = None
        self.gIndices = []
        self.nCycleReacs = 0
        self.r = None
        self.rIndices = []

    @staticmethod
    def computeGamma(S):
        """ perform cycle analysis on matrix S via QR decomposition with column
            pivoting as developed in [1]
        """
        p, _, r = qrpt_decomp2(S)
        n = S.shape[1]     # number of columns in S = number of reactions
        rank = r.shape[0]  # number of rows of R = rank of S
        nRowsG = n-rank
        g = empty((nRowsG, n), S.dtype)
        g[:, :rank] = -r[:rank, rank:].T
        g[:, rank:] = eye(nRowsG)
        return p, r, g

    def getDependent(self):
        """ return list of indices of dependent reactions
        """
        if self.gamma is None:
            return []

        return self.gIndices[-len(self.gamma):]

    def displayGamma(self, reactionNames):
        """ return a list of lines for displaying gamma as linear combinations
        """
        if self.gamma is None:
            return []
        return displayPermutedMatrix(self.gamma, self.gIndices, reactionNames)

    def computeCycles(self, model, internalOnly=True):
        """ perform cycle analysis on the given MetabolicModel

        Results are stored in member variables gamma, gIndices, r, rIndices.
        Afterwards, cycle descriptions can be requested by call of
        getCycleDescriptions().
        """
        if not _pygsl_avail:
            raise ImportError("pygsl not found")

        S = array(model.getStoichiometricMatrix())
        nCols = S.shape[1]
        if internalOnly:
            internalReactions = model.getInternalFluxes()
            submodel = model.getSubmodelByList(internalReactions)
            S_ = array(submodel.getStoichiometricMatrix())
            nColsG = len(submodel)
        else:
            S_ = S
            nColsG = len(model)

        self.gIndices, _, self.gamma = self.computeGamma(S_)

        if internalOnly:
            # Replace indices in submodel with original indices
            for i in range(nColsG):
                self.gIndices[i] = model.reactionDict[submodel.reactions[
                    self.gIndices[i]].name]

        # Round Gamma matrix to 8 digits and truncate near-zero values
        self.gamma = _vmyround(self.gamma)

        # Perform Gauss-Jordan elimination of full matrix (without column
        # pivoting) to get relations between fluxes for all reactions
        # participating in cycles

        # First create permutation of columns of S such that
        # - left-most are the dependent reactions
        # - to the right follow all other reactions participating in any cycle
        # - right-most are the reactions not participating in any cycle
        self.rIndices = empty(nCols, 'i')

        # First store indices of dependent reactions
        nDep = len(self.gamma)
        dependent = self.gIndices[-nDep:]
        self.rIndices[:nDep] = dependent
        # Create mask for column indices already covered by rIndices
        got = zeros(nCols, bool)
        got[dependent] = True

        # Next store indices of independent reactions participating in cycles
        i = nDep
        for j in nonzero(self.gamma[:, :-nDep])[1]:
            index = self.gIndices[j]
            if not got[index]:
                self.rIndices[i] = index
                got[index] = True
                i += 1
        self.nCycleReacs = i

        # Finally add all other reactions
        for index in range(nCols):
            if not got[index]:
                self.rIndices[i] = index
                got[index] = True
                i += 1

        # Permute stoichiometric matrix and perform Gauss-Jordan elimination
        self.r = S[:, self.rIndices]
        rank = gaussJordan(self.r)
        # Remove all-zero rows from R
        self.r = self.r[:rank, :]

        # Also round r matrix to 8 digits and truncate near-zero values
        self.r = _vmyround(self.r)

    def getCycleDescriptions(self, solution, sIndices=None):
        """ generate descriptions of each cycle from Gamma and R, depending on
            the flux distribution given in solution

        Keyword arguments:

        solution -- an FBA solution vector or list of reaction names
        sIndices -- sIndices[i] is the index in the original model corresponding
                    to solution[i] - if not given, assumed as [0..len(solution)]

        Returns: list of FluxCycle objects
        """
        try:
            float(solution[0])
            isString = False
        except ValueError:
            isString = True
        except IndexError:
            return []

        if sIndices == None:
            sIndices = list(range(len(solution)))

        dependent = set(self.getDependent())
        dtype = self.r.dtype
        nColsR = self.r.shape[1]
        cycles = []

        # Build reverse permutations from rIndices and sIndices
        rIndicesRev = empty(nColsR, 'i')
        for i in range(nColsR):
            rIndicesRev[self.rIndices[i]] = i
        sIndicesRev = {}
        for i in range(len(sIndices)):
            sIndicesRev[sIndices[i]] = i

        # Set of all reactions participating in at least one cycle
        cycleIndices = self.rIndices[:self.nCycleReacs]
        allCycleReacs = set(cycleIndices)
#        print "allCycleReacs", allCycleReacs

        # Each row in Gamma corresponds to one cycle
        for rowG in self.gamma:
            cycle = FluxCycle()
            cycle.indices = cycleIndices

            # Only reactions with non-zero coefficient are part of the cycle
            nz = nonzero(rowG)[0]
            nNonzero = len(nz)
            cycle.lhs = empty((nNonzero-1, self.nCycleReacs), dtype)
            if isString:
                cycle.rhs = [""]*(nNonzero-1)
            else:
                cycle.rhs = zeros(nNonzero-1, dtype)
            setCycIndices = set(self.gIndices[nz])

            # Build left-hand side of lin. eq. system from selected columns of R
            # Build right-hand side from masked columns and solution
            i = 0
#            k = 0
#            print "begin Cycle"
#            print " indices", cycle.indices
#            print " processing R"
            for rowR in self.r:
                #                k += 1
                #                print " row %u" % k

                # Only look at non-zero entries of rowR
                nzR = nonzero(rowR)[0]
                rIndicRed = self.rIndices[nzR]
                setRIndicRed = set(rIndicRed)

                # Skip a) rows of R not pertaining to the reactions of the cycle
                #         (i.e. those with all zeros in all selected columns)
                #  and b) rows of R pertaining to more than one cycle
                if (len(setCycIndices & setRIndicRed) == 0 or
                        len((setRIndicRed | setCycIndices) & dependent) > 1):
                    continue

#                print rowR
#                print "  match! indices:", rIndicRed
#                print "  values:", rowR[nzR]
#                sValues = []
#                for j in rIndicRed:
#                    if j in sIndicesRev:
#                        sValues.append(solution[sIndicesRev[j]])
#                    else:
#                        sValues.append(nan)
#                print "  flux values:", sValues
                cycle.lhs[i, :] = rowR[:self.nCycleReacs]
#                print "row in cycle desc.:"
#                print cycle.lhs[i, :]

                if isString:
                    for j in rIndicRed:
                        if j not in allCycleReacs:
                            cycle.rhs[i] += makeCoefNameTerm(cycle.rhs[i],
                                                             -rowR[rIndicesRev[j]], solution[sIndicesRev[j]])
                else:
                    for j in rIndicRed:
                        if j not in allCycleReacs:
                            cycle.rhs[i] -= (rowR[rIndicesRev[j]] *
                                             solution[sIndicesRev[j]])
#                print "rhs:", cycle.rhs[i]

                i += 1
                if i >= nNonzero-1:
                    break

            cycles.append(cycle)

        return cycles
