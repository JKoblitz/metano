"""
This module defines the class LinearProblem, which represents a linear problem
to be solved by an appropriate solver.


This file is part of metano.
Copyright (C) 2010-2014 Alexander Riemer,
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
from openopt import LP
from numpy import dot, inf, isinf
from metano.defines import SolverStatus, glpkToSolverStatus
import pymprog


_GLPK = 'glpk'
_OOGLPK = 'ooglpk'


class SolverError(Exception):
    """ Exception raised if an unknown or inapplicable solver is selected
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class LinearProblem(object):

    # Table of linear problem solvers with corresponding OpenOpt/internal name
    solvers = {'default': 'glpk',
               'lpsolve': 'lpSolve',
               'cvxopt': 'cvxopt_lp',
               'glpk': _GLPK,
               'ooglpk': _OOGLPK}

    def __init__(self, Aeq, beq=None, Aineq=None, bineq=None, lb=None, ub=None,
                 solver=solvers['default']):
        """ initialize the linear problem

        Keyword arguments:

        Aeq         -- matrix of equality constrains (Aeq * v = beq)
        beq         -- right-hand side vector of equality constraints
        Aineq       -- matrix of inequality constrains (Aineq * v <= bineq)
        bineq       -- right-hand side vector of inequality constraints
        lb          -- list of lower bounds (indexed like matrix columns)
        ub          -- list of upper bounds (indexed like matrix columns)
        solver      -- solver to be used for Linear Programming
        """
        """
        with open("parameters.txt", "w") as file:
            string = str(self) + "\n\n Aeq:\n" + str(Aeq) + "\n\n beq:\n" + str(beq) + "\n\n Aineq:\n" + str(Aineq) + "\n\n bineq:\n" + str(bineq) + "\n\n lb:\n" + str(lb) + "\n\n ub:\n" + str(ub) + "\n\n solver:\n" + str(solver) + "\n"
            file.write(string)
"""
        try:
            self.solver = self.solvers[solver.lower()]
        except KeyError:
            raise SolverError("Unknown solver '%s'. Try one of the following:"
                              "\n%s" % (solver, list(self.solvers.keys())))
        self.Aeq = Aeq
        self.Aineq = Aineq
        self.bineq = bineq
        self.istop = None
        self.status = SolverStatus.PREPARED

        # Check matrix and lb/ub vectors
        nRows = len(self.Aeq)
        if not nRows:
            raise ValueError("Error: Matrix Aeq has no rows")
        nCols = len(self.Aeq[0])
        if nCols == 0:
            raise ValueError("Error: Matrix Aeq has no columns")

        if lb is None:
            self.lb = [-inf]*nCols
        elif len(lb) != nCols:
            raise ValueError("Error: Vector of lower bounds has wrong length")
        else:
            self.lb = lb

        if ub is None:
            self.ub = [inf]*nCols
        elif len(ub) != nCols:
            raise ValueError("Error: Vector of upper bounds has wrong length")
        else:
            self.ub = ub

        if beq is not None:
            if len(beq) != nRows:
                raise ValueError("Error: Right-hand side vector of equalities "
                                 "has wrong length")
            self.beq = beq
        else:
            # Default: Right-hand side is zero vector with the same number of
            #          rows as matrix
            self.beq = [0.]*nRows

        if Aineq is None:
            nIneq = 0
        else:
            nIneq = len(Aineq)
            if len(bineq) != nIneq:
                raise ValueError("Error: Matrix and vector dimensions disagree "
                                 "in inequality constraints")
            if nIneq > 0 and len(Aineq[0]) != nCols:
                raise ValueError("Error: Matrix Aineq has wrong number of "
                                 "columns.")

        # Result values
        self.obj_value = 0.  # value of objective function at solution
        self.solution = []  # solution vector

        if self.solver == _GLPK:
            rangeCols = list(range(nCols))

            self.glpkP = pymprog.model('LP')
            self.glpkVar = self.glpkP.var(
                'flux', rangeCols, bounds=(None, None))
            self.glpkP.solver("simplex", msg_lev=pymprog.glpk.GLP_MSG_ERR)

            # Mass balance and other equality constraints
            self.glpkP.st(sum(self.glpkVar[j] * float(Aeq[i][j])
                              for j in rangeCols if Aeq[i][j] != 0.) ==
                          float(beq[i]) for i in range(nRows))

            # Flux bounds constraints
            self.glpkP.st([self.glpkVar[i] >= float(lb[i])
                           for i in rangeCols if not isinf(lb[i])])
            self.glpkP.st([self.glpkVar[i] <= float(ub[i])
                           for i in rangeCols if not isinf(ub[i])])

            # Inequality constraints
            if nIneq:
                self.glpkP.st(sum(self.glpkVar[j] * float(Aineq[i][j])
                                  for j in rangeCols if Aineq[i][j] != 0.)
                              <= float(bineq[i]) for i in range(nIneq))

            # Set to verbose
            # self.glpkP.verb=False

    def setInequalityConstraints(self, Aineq=None, bineq=None):
        """ set inequality constraints for optimization (or unset if None)
        """
        nCols = len(self.lb)
        if Aineq is None:
            nIneq = 0
        else:
            nIneq = len(Aineq)
            if nIneq != len(bineq):
                raise ValueError("Error: Matrix and vector dimensions disagree "
                                 "in inequality constraints")
            if nIneq > 0 and len(Aineq[0]) != nCols:
                raise ValueError("Error: Matrix Aineq has wrong number of "
                                 "columns.")
        self.Aineq = Aineq
        self.bineq = bineq

        if self.solver == _GLPK:
            # Construct new model object (because constraints can't be unset)
            rangeCols = list(range(nCols))

            self.glpkP = pymprog.model('LP')
            self.glpkVar = self.glpkP.var(
                'flux', rangeCols, bounds=(None, None))
            self.glpkP.solver("simplex", msg_lev=pymprog.glpk.GLP_MSG_ERR)
            nRows = len(self.Aeq)

            # Mass balance constraints
            self.glpkP.st(sum(self.glpkVar[j] * float(self.Aeq[i][j])
                              for j in rangeCols if self.Aeq[i][j] != 0.) ==
                          float(self.beq[i]) for i in range(nRows))

            # Flux bounds constraints
            self.glpkP.st([self.glpkVar[i] >= float(self.lb[i])
                           for i in rangeCols if not isinf(self.lb[i])])
            self.glpkP.st([self.glpkVar[i] <= float(self.ub[i])
                           for i in rangeCols if not isinf(self.ub[i])])

            if nIneq:
                # Inequality constraints
                self.glpkP.st(sum(self.glpkVar[j] * float(Aineq[i][j])
                                  for j in rangeCols if Aineq[i][j] != 0.)
                              <= float(bineq[i]) for i in range(nIneq))

    def setObjective(self, obj):
        """ set objective function as coefficient vector

        Keyword arguments:
        obj         -- the linear objective function (coefficients in the linear
                       combination of matrix columns, i.e. fluxes) (negative for
                       maximization, positive for minimization)
        """
        self.obj = obj
        if len(self.obj) != len(self.lb):
            raise ValueError("Error: Objective function vector has wrong "
                             "length")
        if self.solver == _GLPK:
            self.objIndices = [i for i in range(len(obj)) if obj[i] != 0.]
            self.glpkP.min(sum(self.glpkVar[i]*obj[i] for i in self.objIndices),
                           "objective")

    def resetObjective(self):
        self.obj = [0.]*len(self.obj)
        if self.solver == _GLPK:
            self.glpkP.min(sum(self.glpkVar[i]*0. for i in self.objIndices),
                           "objective")

    def minimize(self, **kwargs):
        """ solve the linear problem with the solver given in self.solver

        Returns:
        obj_value, solution

        obj_value -- value of the objective function at the discovered solution
        solution  -- the solution flux vector (indexed like matrix columns)
        """
        if self.solver == _GLPK:
            # Use pyMathProg interface to GLPK solver
            self.glpkP.solve()

            self.istop = self.glpkP.status(str)  # fix for pymprog v1.1
#            print "istop:", self.istop
            self.status = glpkToSolverStatus(self.istop)

            if self.status == SolverStatus.OPTIMAL:
                self.obj_value = self.glpkP.vobj()
                self.solution = [self.glpkVar[i].primal for i in
                                 range(len(self.glpkVar))]
            else:
                self.obj_value = 0.
                self.solution = []
        else:
            # Use OpenOpt interface
            lp = LP(self.obj, A=self.Aineq, Aeq=self.Aeq, b=self.bineq,
                    beq=self.beq, lb=self.lb, ub=self.ub, **kwargs)
            lp.debug = 1
#            lp.iprint=-1 # suppress solver output
            r = lp.solve(self.solver if self.solver != _OOGLPK else 'glpk')

            if r.istop <= 0. or r.ff != r.ff:  # check halting condition
                self.obj_value = 0.
                self.solution = []
                self.istop = r.istop
            else:
                self.obj_value = r.ff
                self.solution = r.xf
                self.istop = r.istop
#            print self.istop

        return self.obj_value, self.solution

    def eval(self, v):
        """ evaluate objective function for vector v
        """
        return dot(self.obj, v)
