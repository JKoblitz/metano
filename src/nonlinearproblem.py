"""
This module defines the class NonLinearProblem, which represents a nonlinear
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
from __future__ import division

from past.utils import old_div
from openopt import NLP
from metano.linearproblem import LinearProblem


class NonLinearProblem(LinearProblem):

    # Table of nonlinear problem solvers with corresponding OpenOpt name
    solvers = {'default': 'algencan',
               'ralg': 'ralg',           # for medium-sized problems;
               #  ill-conditioned
               'cobyla': 'scipy_cobyla',   # recommended
               'algencan': 'algencan',       # recommended
               'slsqp': 'scipy_slsqp',
               'ipopt': 'ipopt',
               'lincher': 'lincher'}        # very primitive

    def __init__(self, Aeq, beq=None, Aineq=None, bineq=None, lb=None, ub=None,
                 solver=solvers['default']):
        """ initialize the nonlinear problem

        Keyword arguments:

        Aeq     -- matrix of equality constrains (Aeq * v = beq)
        beq     -- right-hand side vector of equality constraints
        Aineq   -- matrix of inequality constrains (Aineq * v <= bineq)
        bineq   -- right-hand side vector of inequality constraints
        lb      -- list of lower bounds (indexed like matrix columns)
        ub      -- list of upper bounds (indexed like matrix columns)
        solver  -- solver to be used for Nonlinear Programming
        """
        LinearProblem.__init__(self, Aeq, beq, Aineq, bineq, lb, ub, solver)
        if solver.lower() in LinearProblem.solvers:
            # If problem has been raised from a linear problem, use default
            # solver.
            self.solver = NonLinearProblem.solvers["default"]
        self.obj = None
        self.d_obj = None
        self.nlc = None
        self.d_nlc = None
        self.x0 = [0.]*len(lb)
        self.iterator = None
        self.rlist = None

    def setObjective(self, f):
        self.obj = f

    def resetObjective(self):
        """ shadow LinearProblem.resetObjective """
        pass

    def setObjGrad(self, df):
        self.d_obj = df

    def setNonlinearConstraints(self, nlc):
        self.nlc = nlc

    def setNlcGrad(self, d_nlc):
        self.d_nlc = d_nlc

    def setStartPoint(self, x0):
        if len(x0) != len(self.lb):
            raise ValueError("Error: Starting point has wrong dimension.")
        self.x0 = x0
        self.rlist = None

    def setStartPointIterator(self, iterator):
        self.iterator = iterator

    def minimize(self, **kwargs):
        """ solve the nonlinear problem using OpenOpt

        Returns:
        obj_value, solution

        obj_value -- value of the objective function at the discovered solution
        solution  -- the solution flux vector (indexed like matrix columns)
        """
        if self.iterator is None:
            nlp = NLP(self.obj, self.x0, df=self.d_obj, c=self.nlc,
                      dc=self.d_nlc, A=self.Aineq, Aeq=self.Aeq, b=self.bineq,
                      beq=self.beq, lb=self.lb, ub=self.ub, **kwargs)
            nlp.debug = 1
            nlp.plot = False
            nlp.checkdf()
            if self.nlc is not None:
                nlp.checkdc()

            r = nlp.solve(self.solver)

        else:
            self.rlist = []
            for x0 in self.iterator:
                nlp = NLP(self.obj, x0, df=self.d_obj, c=self.nlc,
                          dc=self.d_nlc, A=self.Aineq, Aeq=self.Aeq,
                          b=self.bineq, beq=self.beq, lb=self.lb, ub=self.ub,
                          **kwargs)
                r = nlp.solve(self.solver)
                if r.istop > 0 and r.ff == r.ff:
                    self.rlist.append(r)

            if self.rlist != []:
                r = min(self.rlist, key=lambda x: x.ff)

        if r.istop <= 0 or r.ff != r.ff:  # check halting condition
            self.obj_value = 0.
            self.solution = []
            self.istop = r.istop
        else:
            self.obj_value = r.ff
            self.solution = r.xf
            self.istop = r.istop

        return self.obj_value, self.solution

    def getHistogram(self, numBins=-1):
        """ compute a histogram over all found solutions (objective function
            values)

            Returns:
            histogram, minval, maxval
        """
        numPoints = len(self.rlist)
        if numPoints == 0:
            return [], 0., 0.

        if numBins == 0:
            raise ValueError("numBins must not be zero")

        if numBins < 0:
            numBins = min(int(round(numPoints / 10.)), 1)

        minval = min(self.rlist, key=lambda x: x.ff).ff
        maxval = max(self.rlist, key=lambda x: x.ff).ff
        if minval == maxval:
            return [numPoints], minval, maxval
        binWidth = (maxval-minval)/float(numBins)

        histogram = [0]*(numBins+1)  # Extra bin for maxval

        for r in self.rlist:
            histogram[int(old_div((r.ff-minval), binWidth))] += 1

        histogram[-2] += histogram[-1]
        del histogram[-1]

        return histogram, minval, maxval

    def eval(self, v):
        """ evaluate objective function for vector v
        """
        return self.d_obj(v)
