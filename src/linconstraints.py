"""
This module defines the class LinearConstraint, which represents a linear
equality or inequality constraint.


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
class LinearConstraint(object):
    """ Represents a single linear equality or inequality constraint

    Member variables:

    lhs  -- list of (coef, name) tuples
    rhs  -- value of right-hand side of equality/inequality (float constant)
    isEq -- if True, interpretation is 'lhs = rhs', else 'lhs < rhs'
    """

    def __init__(self, lhs, sign=None, rhs=None):
        """ create a new linear constraint

        Keyword arguments:

        lhs  -- either a list of (coef, reaction name) tuples or a list of such
                tuples followed by a sign (see below) and a float value -
                explicitly given sign and rhs take precedence over those
                potentially defined in lhs
        sign -- comparison sign (<, >, =, <=, >=, or ==)
        rhs  -- value of right-hand side of equality/inequality (float constant)

        Sign and rhs are optional, defaults are '<' and 0.
        """
        self.clear()
        lhsExtra = []
        for elem in lhs:
            try:
                _, _ = elem
                self.lhs.append(elem)
            except (ValueError, TypeError):
                lhsExtra.append(elem)

        if lhsExtra:
            if not sign:
                sign = lhsExtra[0]
            if len(lhsExtra) > 1:
                if not rhs:
                    rhs = lhsExtra[1]

        if rhs:
            self.rhs = rhs

        # Parse sign
        if sign:
            sign = sign[0]
            if sign == '=':
                self.isEq = True

            elif sign == '>':
                # Flip inequality by negating both signs
                self.rhs = -self.rhs
                self.lhs = [(-coef, name) for coef, name in self.lhs]


    def clear(self):
        self.lhs = []
        self.rhs = 0.
        self.isEq = False


    def getSign(self):
        return '=' if self.isEq else '<'


    def __str__(self):
        s = ""
        for coef, name in self.lhs:
            if not s:
                s = "%g %s" % (coef, name)
            else:
                if coef < 0.:
                    s+= " - %g %s" % (-coef, name)
                else:
                    s+= " + %g %s" % (coef, name)
        return "%s %s %g" % (s, self.getSign(), self.rhs)


    def __repr__(self):
        return "LinearConstraint(%r, %r, %r)" % (self.lhs, self.getSign(),
                                                 self.rhs)
