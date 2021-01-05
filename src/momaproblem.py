"""
This module defines the class MomaProblem, which represents a quadratic problem
to be solved in MOMA analysis (Maximization of Metabolic Adjustment).

The MomaProblem is a quadratic optimization problem with objective function
    Z = (v-v_wt)'*(v-v_wt) = v'*v - 2*v_wt'*v + v_wt'*v_wt,
i.e. the goal is to minimize the quadratic distance to the FBA solution for the
wildtype.

OpenOpt demands a formulation in terms of a matrix H and a vector f:
    Z_oo = 0.5 v'*H*v + f'*v

For H = Identity and f = -v_wt, minimization of Z_oo is equivalent to
minimization of Z, because the constant v_wt'*v_wt can be omitted and the
function can be scaled by a factor.

If both CVXMOD and CVXPY are installed, the CVXMOD interface will be used
preferably, as it is faster. Only if CVXMOD fails to converge, CVXPY is used.


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

from builtins import range
from numpy import array, eye, inf, isnan, diag, sqrt, linalg
from metano.defines import SolverStatus, cvxToSolverStatus, cvxpyToSolverStatus
from metano.linearproblem import SolverError
from metano.quadraticproblem import QuadraticProblem, _cvxmod_avail, _cvxpy_avail
if _cvxmod_avail:
    import cvxmod
if _cvxpy_avail:
    import cvxpy


class MomaProblem(QuadraticProblem):

    matrix_tol = 1E-5     # tolerance for matrix constraints
    bounds_tol = 1E-6     # tolerance for bounds constraints
    BIG_NUMBER = 100000

    def __init__(self, Aeq, beq=None, Aineq=None, bineq=None, lb=None, ub=None,
                 solver=QuadraticProblem.solvers['default'], solution=None,
                 weights=None):
        """ initialize the quadratic problem

        Keyword arguments:

        Aeq         -- matrix of equality constrains (Aeq * v = beq)
        beq         -- right-hand side vector of equality constraints
        Aineq       -- matrix of inequality constrains (Aineq * v <= bineq)
        bineq       -- right-hand side vector of inequality constraints
        lb          -- list of lower bounds (indexed like matrix columns)
        ub          -- list of upper bounds (indexed like matrix columns)
        solver      -- QP or NLP solver to be used for optimization
        solution    -- the solution for the wildtype
        weights     -- weight vector
        """
        print("Momaproblem", solver)
        QuadraticProblem.__init__(self, Aeq, beq, Aineq, bineq, lb, ub, solver)
        if weights is None:
            self.weights = None
            self.obj.H = eye(len(lb))
        else:
            self.weights = array(weights)
            if len(weights) != len(lb):
                raise ValueError("Error: Weight vector has wrong length")
            self.obj.H = diag(self.weights)

        if solution is not None:
            self.obj.f = -array(solution)

    def setWtSolution(self, solution):
        """ set the solution for the wildtype against which to perform MOMA in
            the objective function definition of the underlying QuadraticProblem
        """
        self.obj.f = -array(solution)

    def minimize(self, **kwargs):
        if self.solver.startswith("cvx"):
            if not _cvxmod_avail and not _cvxpy_avail:
                raise SolverError("Solver '%s' not available. See installation "
                                  "instructions on http://metano.tu-bs.de" %
                                  self.solver)
            nCols = len(self.lb)
            lb = [-self.BIG_NUMBER]*nCols
            ub = [self.BIG_NUMBER]*nCols
            for i in range(nCols):
                if self.lb[i] > -inf:
                    lb[i] = self.lb[i]-self.bounds_tol
                if self.ub[i] < inf:
                    ub[i] = self.ub[i]+self.bounds_tol

            doCvxPy = self.solver == "cvxpy_only" or not _cvxmod_avail
            if _cvxmod_avail and (self.solver == "cvxopt_qp" or
                                  self.solver == "cvxmod_only"):
                self._cvxmod_minimize(lb, ub, **kwargs)
                if self.status == SolverStatus.UNKNOWN:
                    doCvxPy = True

            if doCvxPy and _cvxpy_avail:
                self._cvxpy_minimize(lb, ub, **kwargs)

            return self.obj_value, self.solution
        else:  # use OpenOpt solver
            if self.weights is not None:
                self.obj.f *= self.weights
            solution = QuadraticProblem.minimize(self, **kwargs)
            if self.weights is not None:
                self.obj.f /= self.weights
            return solution

    def _cvxmod_minimize(self, lb, ub, **kwargs):
        """ solve quadratic problem using CVXMOD interface

        Keyword arguments:

        lb, ub  -- vectors of lower and upper bounds (potentially modified by
                   added tolerance)

        Modified member variables:

        status, solution, obj_value

        Returns: nothing
        """
        # shorter names
        v = self.cvxmodV
        A = self.cvxmodMatrix
        minus_w = cvxmod.matrix(self.obj.f)  # negative wildtype solution
        lb = cvxmod.matrix(lb)
        ub = cvxmod.matrix(ub)
        if self.weights is not None:
            weights = cvxmod.matrix(diag(sqrt(self.weights)))
        if self.Aineq is not None:
            Aineq = cvxmod.matrix(array(self.Aineq))
            bineq = cvxmod.matrix(self.bineq)
            if self.weights is None:
                p = cvxmod.problem(cvxmod.minimize(cvxmod.norm2(v+minus_w)),
                                   [cvxmod.abs(A*v) < self.matrix_tol,
                                    Aineq*v <= bineq+self.matrix_tol,
                                    v >= lb, v <= ub])
            else:
                p = cvxmod.problem(cvxmod.minimize(cvxmod.norm2(weights *
                                                                (v+minus_w))), [cvxmod.abs(A*v) < self.matrix_tol,
                                                                                Aineq*v <= bineq+self.matrix_tol,
                                                                                v >= lb, v <= ub])
        else:
            if self.weights is None:
                p = cvxmod.problem(cvxmod.minimize(cvxmod.norm2(v+minus_w)),
                                   [cvxmod.abs(A*v) < self.matrix_tol,
                                    v >= lb, v <= ub])
            else:
                p = cvxmod.problem(cvxmod.minimize(cvxmod.norm2(weights *
                                                                (v+minus_w))), [cvxmod.abs(A*v) < self.matrix_tol,
                                                                                v >= lb, v <= ub])

        self.status = cvxToSolverStatus(p.solve())

        if not v.value:
            self.solution = []
        else:
            self.solution = array(list(v.value))
        try:
            self.obj_value = p.value
        except cvxmod.OptvarValueError:
            self.obj_value = inf

    def _cvxpy_minimize(self, lb, ub, **kwargs):
        """ solve quadratic problem using CVXPY interface

        Keyword arguments:

        lb, ub  -- vectors of lower and upper bounds (potentially modified by
                   added tolerance)

        Modified member variables:

        status, solution, obj_value

        Returns: nothing
        """
        # shorter names
        v = self.cvxpyV
        A = self.cvxpyMatrix
        minus_w = array(self.obj.f).T  # negative wildtype solution
        lb = array(lb).T
        ub = array(ub).T
        if self.weights is not None:
            weights = array(diag(sqrt(self.weights)))
        if self.Aineq is not None:
            Aineq = array(self.Aineq)
            bineq = array(self.bineq)
            p = cvxpy.Problem(cvxpy.Minimize(cvxpy.norm2(v+minus_w)),
                              [cvxpy.abs(A*v) <= self.matrix_tol,
                               (Aineq*v, bineq) <= v, v >= lb,
                               v <= ub])
        else:
            p = cvxpy.Problem(cvxpy.Minimize(cvxpy.norm2(v+minus_w)),
                              [cvxpy.abs(A*v) <= self.matrix_tol,
                               v >= lb, v <= ub])

        # Configure program ('abstol', 'feastol', 'reltol', 'maxiters')
        """for k in kwargs:
            try:
                p.options[k] = kwargs[k]
            except KeyError:
                continue"""
        self.obj_value = p.solve(solver="CVXOPT", reltol=kwargs['reltol'],
                                 abstol=kwargs['abstol'], feastol=kwargs['feastol'])
        self.status = cvxpyToSolverStatus(self.obj_value)

        if len(v.value) == 0 or isnan(v.value[0]):
            self.solution = []
        else:
            self.solution = array(v.value.T)[:]

    def eval(self, v):
        """ evaluate objective function for vector v
        """
        if self.weights is None:
            return linalg.norm(array(v)+self.obj.f)
        else:
            return linalg.norm(diag(sqrt(self.weights))*(array(v)+self.obj.f))
