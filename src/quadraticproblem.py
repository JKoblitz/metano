"""
This module defines the class QuadraticProblem, which represents a quadratic
problem to be solved by an appropriate solver.


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

from builtins import object
from openopt import QP
from numpy import array
from metano.linearproblem import LinearProblem

try:
    from cvxmod import matrix as cvxmod_matrix, optvar
    _cvxmod_avail = True
except ImportError:
    _cvxmod_avail = False
try:
    from cvxpy import *
    _cvxpy_avail = True
except ImportError:
    _cvxpy_avail = False


class QuadraticProblem(LinearProblem):
    # Table of quadratic problem solvers with corresponding OpenOpt name
    solvers = {'default': 'cvxopt_qp',          # cvxopt = try cvxmod first
               'cvxopt': 'cvxopt_qp',  # (fast) - if no convergence,
               'cvxopt_qp': 'cvxopt_qp',  # try cvxpy (slower)
               'cvxmod': 'cvxmod_only',
               'cvxpy': 'cvxpy_only',
               # Or use an NLP solver:
               'ralg': 'nlp:ralg',           # for medium-sized problems;
               #   ill-conditioned
               'cobyla': 'nlp:scipy_cobyla',   # recommended (for NLP)
               'algencan': 'nlp:algencan',       # recommended (for NLP)
               'slsqp': 'nlp:scipy_slsqp',
               'ipopt': 'nlp:ipopt',
               'lincher': 'nlp:lincher'}        # very primitive

    @staticmethod
    def isCvxQpSolver(solver):
        return solver.startswith("cvx") or solver == "default"

    class OoQp(object):
        """ class for objective function in OpenOpt's QP format:

        Z(v) = 0.5*v'*H*v + f*v with

            H - matrix of coefficients of quadratic terms and
            f - vector of coefficients of linear terms
        """

        def __init__(self, H=None, f=None):
            self.H = H
            self.f = f

    def __init__(self, Aeq, beq=None, Aineq=None, bineq=None, lb=None, ub=None,
                 solver=solvers['default']):
        """ initialize the quadratic problem

        Keyword arguments:

        Aeq     -- matrix of equality constrains (Aeq * v = beq)
        beq     -- right-hand side vector of equality constraints
        Aineq   -- matrix of inequality constrains (Aineq * v <= bineq)
        bineq   -- right-hand side vector of inequality constraints
        lb      -- list of lower bounds (indexed like matrix columns)
        ub      -- list of upper bounds (indexed like matrix columns)
        solver  -- solver to be used for Quadratic Programming
        """
        LinearProblem.__init__(self, Aeq, beq, Aineq, bineq, lb, ub, solver)

        self.obj = QuadraticProblem.OoQp()
        if _cvxmod_avail:
            self.cvxmodMatrix = cvxmod_matrix(array(Aeq))
            self.cvxmodV = optvar('v', self.cvxmodMatrix.size[1], 1)
        if _cvxpy_avail:
            self.cvxpyMatrix = array(Aeq)
            self.cvxpyV = Variable(self.cvxpyMatrix.shape[1], 'v', 1)

    def setObjective(self, quad_coef, lin_coef):
        """ set the objective function in OpenOpt's QP format:

        Z(v) = 0.5*v'*H*v + f*v with

            H = quad_coef, matrix of coefficients of quadratic terms and
            f = lin_coef, vector of coefficients of linear terms
        """
        self.obj = QuadraticProblem.OoQp(quad_coef, lin_coef)

    def resetObjective(self):
        """ shadow LinearProblem.resetObjective """
        pass

    def minimize(self, **kwargs):
        """ solve the quadratic problem using OpenOpt

        Returns:
        obj_value, solution

        obj_value -- value of the objective function at the discovered solution
        solution  -- the solution flux vector (indexed like matrix columns)
        """
        qp = QP(self.obj.H, self.obj.f, A=self.Aineq, Aeq=self.Aeq,
                b=self.bineq, beq=self.beq, lb=self.lb, ub=self.ub, **kwargs)
        qp.debug = 1
        r = qp.solve(self.solver)

        if r.istop <= 0 or r.ff != r.ff:  # check halting condition
            self.obj_value = 0
            self.solution = []
            self.istop = r.istop
        else:
            self.obj_value = r.ff
            self.solution = r.xf
            self.istop = r.istop

        return self.obj_value, self.solution
