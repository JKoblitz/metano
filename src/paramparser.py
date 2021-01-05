"""
This module defines the class ParamParser, a parser for scenario files for FBA.


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

# TODO Release:
# - remove SAMESIGN/DIFFSIGN keywords & getNLConstraints()

from builtins import range
from past.utils import old_div
from builtins import object
from metano.defines import Verbosity
from metano.linconstraints import LinearConstraint
from numpy import array, inf, any
import re


class LinExprParser(object):
    """ Parser for linear expressions of the form
        <coefficient> <name> [+|- <coefficient> <name>]*

    This is a static class.
    Note: As +/- and numbers can be parts of names in metano, spaces around
          operators and numbers are important!
    """

    _operators = ("+-*")
    _escOp = {"+": "\+", "*": "\*"}  # dict for escaping operators for reg-exes
    _floatRegexBase = (r"(?:^|\s)"  # space or beginning of string
                       r"([-+]?"     # sign (optional)
                       # arbitrary number of digits (optional)
                       # either (digit [dot]) or (dot digit)
                       # arbitrary number of digits (optional)
                       r"(?:[0-9]+\.?[0-9]*|[0-9]*\.[0-9]+)"
                       r"(?:[eE][-+]?[0-9]+)?)")  # exponent (e/E [+/-] digits)
    _floatRegexSpc = _floatRegexBase + r"(?:$|\s)"  # space or end of string
    _floatSpcPat = re.compile(_floatRegexSpc)
    # regex matching a colon-marked float
    _floatMarked = _floatRegexBase+":"
    _term = "(?:"+_floatMarked+r"(?:\s*\*:)?\s*)?([^\s:]+)(?:$|\s+)"
    _term1 = re.compile(_term)
    _term2ff = re.compile(r"([+-]):\s*"+_term)

    @staticmethod
    def parseString(s):
        sOrig = s

        # Preprocessing: Mark all numbers and mathematical operators that are
        # not part of words with the illegal character ':'
        for op in LinExprParser._operators:
            s = (" "+op+": ").join(re.split("\s"+LinExprParser._escOp.get(op,
                                                                          op)+"(?:\s|$)", s))
        floatHits = list(LinExprParser._floatSpcPat.finditer(s))
        floatPos = [i.end() for i in floatHits[::-1]]
        for i in floatPos:
            s = s[:i].rstrip()+": "+s[i:]

        # Parsing: Compile list of (sign, coef, name) tuples
        t = []
        m = LinExprParser._term1.match(s)
        if not m:
            raise ValueError("string '%s' is not a linear expression" % sOrig)
        t.append(("+",)+m.groups())
        s = s[m.end():]

        while s:
            m = LinExprParser._term2ff.match(s)
            if not m:
                raise ValueError("string '%s' is not a linear expression" %
                                 sOrig)
            t.append(m.groups())
            s = s[m.end():]

        return LinExprParser._evalGroups(t)

    @staticmethod
    def _evalGroups(listOfTuples):
        """ decode a list of linear terms given as a (sign, number, name) tuple,
        where sign is '-' or '+', number is a string or None, and name is a
        string; 'number' may itself carry a sign

        Returns: list of (coefficient, name) tuples with float coefficients
        """
        tuples = []
        for sign, number, name in listOfTuples:
            if not number:
                number = 1.

            if sign == '-':
                tuples.append((-float(number), name))
            else:
                tuples.append((float(number), name))

        return tuples


class LinIneqParser(object):
    """ Parser for linear equations and inequalities of the form
        <linear expression> <|>|=|<=|>=|== <constant>

    This is a static class.
    Note: As +/- and numbers can be parts of names in metano, spaces around
          operators and numbers are important!
    """
    _rel = "(=|<|>)=?"

    @staticmethod
    def parseString(s):
        parts = re.split(LinIneqParser._rel, s)
        if len(parts) != 3:
            raise ValueError("string '%s' is not a linear inequality")
        return LinearConstraint(LinExprParser.parseString(parts[0]) +
                                [parts[1], float(parts[2])])


class ParamParser(object):
    """Parser for scenario files for flux balance analysis.

    The scenario file defines the constraints and objective function and thus
    the optimization problem to be solved. The following lines are recognized:

    OBJ MAX|MIN <linear expression> | PER_FLUX(<linear expression>)
        -- defines the objective function, either as a linear combination of
           reactions (i.e. fluxes through these reactions) or by a built-in
           function in such a linear expression; MAX|MIN specifies the direction
           of optimization (minimization or maximization)
    SOLVER         -- name of preferred LP/QP/NLP solver to be used in FBA
    NUM_ITER <value>           -- number of iterations (in case of NLP)
    LB <reaction-name> <value> -- lower bound for flux through reaction
    UB <reaction-name> <value> -- upper bound for flux through reaction
    SAMESIGN <reaction-name-1> <reaction-name-2>
        -- NL constraint: fluxes through reactions 1&2 must have same sign
    DIFFSIGN <reaction-name-1> <reaction-name-2>
        -- NL constraint: fluxes through reactions 1&2 must have opposite signs

    Moreover, linear inequality constraints can be written as equations or
    inequalities, e.g.: R001 + 2 R002 < 10

    Keywords and MAX|MIN are case-insensitive.

    Available solvers (names are case-insensitive):
        LP: lpsolve (default), gplk, cvxopt (doesn't seem to work)
    """
    ALL_SYMB = "$ALL"

    obj_hint = "OBJ MAX|MIN <linear expression>|PER_FLUX(<linear expression>)"
    solver_hint = "SOLVER <name>"
    numiter_hint = "NUM_ITER <value>"
    lb_hint = "LB <reaction-name> <value>"
    ub_hint = "UB <reaction-name> <value>"
    samesign_hint = "SAMESIGN <reaction-name-1> <reaction-name-2>"
    diffsign_hint = "DIFFSIGN <reaction-name-1> <reaction-name-2>"

    # Keywords that must be present (at least once)
    required_keywords = ["OBJ"]
    # Keywords that must be unique, i.e. may be present at most once
    unique_keywords = ["OBJ", "SOLVER"]

    def __init__(self, commentSign='#'):
        """ initializes the ParamParser

        Keyword arguments:

        commentSign -- comment sign (if a line contains this sign, the sign and
                       everything following it are ignored)
        """
        self.clear()
        self.comment_sign = commentSign

    def clear(self):
        """ clears all internal variables, warnings and info messages
        """
        self.maxmin = False  # True: maximize, False: minimize objective function
        self.obj_name = ""
        self.solver = ""
        self.numIter = -1  # Negative number -> not set
        self.lb = {}  # dictionary {reaction_name : value}
        self.ub = {}  # - '' -
        self.samesign_pairs = set()
        self.diffsign_pairs = set()
        self.lin_constraints = []  # list of LinearConstraint objects
        # list of (line number, original line) pairs for linear inequalities
        self.lin_constraints_orig = []

        # Flags for unique (i.e. may occur only once) and required entries
        # - dictionary {keyword : got-flag}
        self.got = dict([(kw, False) for kw in self.unique_keywords])
        for kw in self.required_keywords:
            if kw not in self.got:
                self.got[kw] = False

        # List of warnings and info messages generated during parsing,
        # elements are pairs (level, msg)
        self.messages = []

    def getMessages(self, level=Verbosity.INFO):
        """ returns a list of all messages at or above the given level of
            severity
            Levels are defined in defines.py.
        """
        return [x[1] for x in self.messages if x[0] <= level]

    def parse(self, filename):
        """ parse the scenario file identified by filename
            -- wrapper for parseByHandle()
        """
        with open(filename) as f:
            return self.parseByHandle(f)

    def parseByHandle(self, f):
        """ parse the scenario file given as a file object (must be open)

        Returns:
        maxmin, obj_name, solver, numIter, spec_lb, spec_ub

        maxmin      -- Goal is to maximize (True) or minimize(False) obj. func.
        obj_name    -- name (or equation) of objective function
        solver      -- name of preferred LP/QP/NLP solver (may be empty)
        numIter     -- number of iterations in case of nonlinear optimization
        spec_lb     -- dictionary of specific lower bounds (by reaction name)
        spec_ub     -- dictionary of specific upper bounds (- '' -)
        """
        # Table of functions for parsing lines by keyword (1st word of line)
        line_parsers = {"OBJ": self.parseObj,
                        "SOLVER": self.parseSolver,
                        "NUM_ITER": self.parseNumIter,
                        "LB": self.parseLb,
                        "UB": self.parseUb,
                        "SAMESIGN": self.parseSameSign,
                        "DIFFSIGN": self.parseDiffSign}

        line_no = 0
        for line in f:
            line_no += 1
            comment_pos = line.find(self.comment_sign)
            line = line[:comment_pos].strip(
            ) if comment_pos >= 0 else line.strip()
            if line == "":
                continue  # Skip blank lines

            # Try to interpret line as equation/inequality definition
            try:
                self.lin_constraints.append(LinIneqParser.parseString(line))
                self.lin_constraints_orig.append((line_no, line))
                isInequality = True
            except ValueError:
                isInequality = False

            if not isInequality:
                try:
                    [kw_original, rest] = line.split(None, 1)
                except ValueError:
                    kw_original = line
                    rest = ""

                keyword = kw_original.upper()
                if keyword in self.unique_keywords:
                    self.checkUnique(keyword, line_no)
                try:
                    line_parsers[keyword](rest, line_no)
                except KeyError:
                    msg = self.error(line_no, "Illegal keyword '%s' or "
                                     "unable to parse linear expression "
                                     "'%s'" % (kw_original, line))
                    raise SyntaxError(msg)

        # Check for presence of required definitions
        for kw in self.required_keywords:
            if not self.got[kw]:
                msg = self.error(0, "Missing %s definition" % kw)
                raise SyntaxError(msg)

        # Simplify linear constraints, convert simple ones to LB/UB definitions
        self._simplifyLinIneqs()
        # Make sure that LB <= UB for all reactions with both LB and UB
        self.checkBounds()

        return (self.maxmin, self.obj_name, self.solver, self.numIter, self.lb,
                self.ub)

    @staticmethod
    def linExprToVector(linExpr, reactions, factor=1., nCols=-1, strict=True):
        """ convert the given linear expression to a coefficient vector

        Keyword arguments:

        linExpr     -- a) list of (coefficient, name) tuples or
                       b) a single reaction name
        reactions   -- a) dictionary of all reactions { name : matrix column }
                       or b) dictionary of split reactions { name : tuple }
        factor      -- factor with which to multiply the result vector
        nCols       -- number of matrix columns or -1 for len(reactions)
        strict      -- if True raise KeyError if a reaction is not found in
                       dict; if False treat as zero (skip)

        Returns: list of coefficients, indexed like reactions dictionary
        """
        if nCols < 0:
            nCols = len(reactions)
        if nCols == 0:
            return []  # nothing to do
        coefVec = [0.]*nCols

        try:
            iter(list(reactions.values())[0])
            isSplit = True
        except TypeError:
            isSplit = False
        if isSplit:
            # Case 1: Flux variables are split into non-negative components
            if isinstance(linExpr, str):
                # 1.a) linExpr is a single reaction name
                try:
                    indexPos, indexNeg = reactions[linExpr]
                    if indexPos is not None:
                        coefVec[indexPos] = factor
                    if indexNeg is not None:
                        coefVec[indexNeg] = -factor
                except KeyError:
                    if strict:
                        raise KeyError("Reaction %s not found in reactions dict"
                                       % linExpr)
                return coefVec

            # 1.b) true linear expression
            for coef, name in linExpr:
                if name.upper() == ParamParser.ALL_SYMB:
                    # $ALL entry defines a global base coefficient (usually 1)
                    for indexPos, indexNeg in reactions:
                        if indexPos is not None:
                            coefVec[indexPos] += coef*factor
                        if indexNeg is not None:
                            coefVec[indexNeg] -= coef*factor
                else:
                    try:
                        indexPos, indexNeg = reactions[name]
                        if indexPos is not None:
                            coefVec[indexPos] += coef*factor
                        if indexNeg is not None:
                            coefVec[indexNeg] -= coef*factor
                    except KeyError:
                        if strict:
                            raise KeyError("Reaction %s not found in reactions "
                                           "dict" % name)

        else:
            # Case 2: Flux variables are not split, i.e. each matrix column
            #         corresponds to one reaction
            if isinstance(linExpr, str):
                # 2.a) linExpr is a single reaction name
                try:
                    coefVec[reactions[linExpr]] = factor
                except KeyError:
                    if strict:
                        raise KeyError("Reaction %s not found in reactions dict"
                                       % linExpr)
                return coefVec

            # 2.b) true linear expression
            for coef, name in linExpr:
                if name.upper() == ParamParser.ALL_SYMB:
                    # $ALL entry defines a global base coefficient (usually 1)
                    for i in range(nCols):
                        coefVec[i] += coef*factor
                else:
                    try:
                        coefVec[reactions[name]] += coef*factor
                    except KeyError:
                        if strict:
                            raise KeyError("Reaction %s not found in reactions "
                                           "dict" % name)
        return coefVec

    @staticmethod
    def linConstraintsToVectors(constraints, reactions, nCols=-1):
        """ get the linear constraints as (coefficient vector, right-hand side)
            pairs

        Any reaction not found in 'reactions' is skipped (treated as zero).

        Keyword arguments:

        constraints -- list of LinearConstraint objects
        reactions   -- a) dictionary of all reactions { name : matrix column }
                       or b) dictionary of split reactions { name : tuple }
        nCols       -- number of matrix columns or -1 for len(reactions)

        Returns: (eqs, ineqs) with

        eqs       -- list of (coefficient vector, right-hand side) pairs for
                     the equality constraints
        ineqs     -- list of - '' - for the inequality constraints
        """
        eqs = []   # Equality constraints
        ineqs = []  # Inequality constraints
        for elem in constraints:
            coefVec = ParamParser.linExprToVector(elem.lhs, reactions,
                                                  nCols=nCols, strict=False)
            if elem.isEq:
                if any(coefVec):
                    eqs.append((coefVec, elem.rhs))
                elif abs(elem.rhs) > 1e-6:
                    raise ValueError("constraint '%s' is always false" % elem)
            else:
                if any(coefVec):
                    ineqs.append((coefVec, elem.rhs))
                elif elem.rhs < 0.:
                    raise ValueError("constraint '%s' is always false" % elem)

        return eqs, ineqs

    @staticmethod
    def convertObjFuncToLinVec(objStr, reactions, nCols=-1, maxmin=False):
        """ construct the coefficient vector for the given linear objective
            function

        Keyword arguments:

        objStr     -- objective function definition (given as string)
        reactions  -- dictionary of all reactions { name : matrix column } or
                      { name : tuple of indices } for split fluxes
        nCols      -- number of matrix columns (or -1 for len(reactions))
        maxmin     -- True for maximization, False for minimization

        Return:  vector of coefficients (same size as flux vector)
        """
        maxmin_factor = -1. if maxmin == True else 1.
        if objStr in reactions:
            objExpr = objStr
        else:
            try:
                objExpr = LinExprParser.parseString(objStr)
            except ValueError:
                raise ValueError("Objective function '%s' is not a linear "
                                 "expression" % objStr)
        return ParamParser.linExprToVector(objExpr, reactions, maxmin_factor,
                                           nCols)

    def getNLConstraints(self, reactions):
        """ build lists of nonlinear constraints (lambda functions) from
            samesign_pairs and diffsign_pairs

        Keyword arguments:

        reactions   -- dictionary of all reactions { name : matrix column }

        Returns:
        nlc         -- list of nonlinear constraints (lambda functions)
        nlc_grad    -- gradients of nonlinear constraints
        """
        if len(self.samesign_pairs) == len(self.diffsign_pairs) == 0:
            return [], []

        numReactions = len(reactions)
        nlc, nlc_grad = [], []
        nonexisting = set()
        for r1, r2 in self.samesign_pairs:
            skip1 = r1 not in reactions
            skip2 = r2 not in reactions
            if skip1:
                if r1 not in nonexisting:
                    msg = ("Warning: SAMESIGN defined for non-existing reaction"
                           " %s" % r1)
                    self.messages.append((Verbosity.WARNING, msg))
                    nonexisting.add(r1)
            if skip2:
                if r2 not in nonexisting:
                    msg = ("Warning: SAMESIGN defined for non-existing reaction"
                           " %s" % r2)
                    self.messages.append((Verbosity.WARNING, msg))
                    nonexisting.add(r2)

            if skip1 or skip2:
                continue

            r1index, r2index = reactions[r1], reactions[r2]
            nlc.append(lambda v: -v[r1index]*v[r2index])

            if r1index > r2index:
                # Make sure that r1index < r2index
                tmp = r1index
                r1index = r2index
                r2index = tmp

            nlc_grad.append(lambda v: array(
                [0. for _ in range(r1index)]+[-v[r2index]] +
                [0. for _ in range(r1index+1, r2index)]+[-v[r1index]] +
                [0. for _ in range(r2index+1, numReactions)]))

        nonexisting.clear()
        for r1, r2 in self.diffsign_pairs:
            skip1 = r1 not in reactions
            skip2 = r2 not in reactions
            if skip1:
                if r1 not in nonexisting:
                    msg = ("Warning: DIFFSIGN defined for non-existing reaction"
                           " %s" % r1)
                    self.messages.append((Verbosity.WARNING, msg))
                    nonexisting.add(r1)
            if skip2:
                if r2 not in nonexisting:
                    msg = ("Warning: DIFFSIGN defined for non-existing reaction"
                           " %s" % r2)
                    self.messages.append((Verbosity.WARNING, msg))
                    nonexisting.add(r2)

            if skip1 or skip2:
                continue

            r1index, r2index = reactions[r1], reactions[r2]
            nlc.append(lambda v: v[r1index]*v[r2index])

            if r1index > r2index:
                # Make sure that r1index < r2index
                tmp = r1index
                r1index = r2index
                r2index = tmp

            nlc_grad.append(lambda v: array(
                [0. for _ in range(r1index)]+[v[r2index]] +
                [0. for _ in range(r1index+1, r2index)]+[v[r1index]] +
                [0. for _ in range(r2index+1, numReactions)]))
        return nlc, nlc_grad

    def parseObj(self, line, line_no=0):
        """ parses the OBJ line, which defines the objective function
        """
        try:
            [mm, name] = line.split(None, 1)
        except ValueError:
            msg = self.error(
                line_no, "Malformed OBJ definition", self.obj_hint)
            raise SyntaxError(msg)

        try:
            self.maxmin = {"MAX": True, "MIN": False}[mm.upper()]
        except KeyError:
            msg = self.error(
                line_no, "Malformed OBJ definition", self.obj_hint)
            raise SyntaxError(msg)

        self.obj_name = name

    def parseSolver(self, line, line_no=0):
        """ parses the SOLVER line, which defines the preferred LP/QP/NLP solver
        """
        self.solver = line.strip()
        if self.solver == "":
            msg = self.error(line_no, "Malformed SOLVER definition",
                             self.solver_hint)
            raise SyntaxError(msg)

    def parseNumIter(self, line, line_no=0):
        """ parses a NUM_ITER line, which defines the number of iterations in
            case of nonlinear optimization
        """
        if line == "" or len(line.split(None, 1)) > 1:
            msg = self.error(line_no, "Malformed NUM_ITER definition",
                             self.numiter_hint)
            raise SyntaxError(msg)

        try:
            self.numIter = int(line)
        except ValueError:
            msg = self.error(line_no, "Illegal integer value in NUM_ITER "
                             "definition", self.numiter_hint)
            raise SyntaxError(msg)

    def parseLb(self, line, line_no=0):
        """ parses an LB line, which defines a specific lower bound
        """
        try:
            [name, value] = line.rsplit(None, 1)
        except ValueError:
            msg = self.error(line_no, "Malformed LB definition", self.lb_hint)
            raise SyntaxError(msg)

        if name in self.lb:
            msg = self.error(line_no, "Duplicate LB definition for %s "
                             "encountered" % name)
            raise SyntaxError(msg)

        try:
            self.lb[name] = float(value)
        except ValueError:
            msg = self.error(line_no, "Illegal floating-point value in LB "
                             "definition", self.lb_hint)
            raise SyntaxError(msg)

    def parseUb(self, line, line_no=0):
        """ parses a UB line, which defines a specific upper bound
        """
        try:
            [name, value] = line.rsplit(None, 1)
        except ValueError:
            msg = self.error(line_no, "Malformed UB definition", self.ub_hint)
            raise SyntaxError(msg)

        if name in self.ub:
            msg = self.error(line_no, "Duplicate UB definition for %s "
                             "encountered" % name)
            raise SyntaxError(msg)

        try:
            self.ub[name] = float(value)
        except ValueError:
            msg = self.error(line_no, "Illegal floating-point value in UB "
                             "definition", self.ub_hint)
            raise SyntaxError(msg)

    def parseSameSign(self, line, line_no=0):
        """ parses a SAMESIGN line, which defines a nonlinear constraint that
            the fluxes through two reactions must have the same sign
        """
        try:
            r1, r2 = line.split(None, 1)
        except ValueError:
            msg = self.error(line_no, "Malformed SAMESIGN definition",
                             self.samesign_hint)
            raise SyntaxError(msg)

        r2 = r2.rstrip()
        if (r1, r2) in self.diffsign_pairs or (r2, r1) in self.diffsign_pairs:
            msg = self.error(line_no, "SAMESIGN definition for reactions "
                             "previously defined as DIFFSIGN. Error")
            raise SyntaxError(msg)

        addToSet = True
        if (r1, r2) in self.samesign_pairs or (r2, r1) in self.samesign_pairs:
            msg = "Warning: Duplicate SAMESIGN definition in line %u" % line_no
            self.messages.append((Verbosity.WARNING, msg))
            addToSet = False

        if r1 == r2:
            msg = ("Warning: SAMESIGN definition in line %u ignored because "
                   "reactions are the same" % line_no)
            self.messages.append((Verbosity.WARNING, msg))
            addToSet = False

        if addToSet:
            self.samesign_pairs.add((r1, r2))

    def parseDiffSign(self, line, line_no=0):
        """ parses a DIFFSIGN line, which defines a nonlinear constraint that
            the fluxes through two reactions must have opposite signs
        """
        try:
            r1, r2 = line.split(None, 1)
        except ValueError:
            msg = self.error(line_no, "Malformed DIFFSIGN definition",
                             self.diffsign_hint)
            raise SyntaxError(msg)

        r2 = r2.rstrip()
        if r1 == r2:
            msg = self.error(line_no, "Reactions in DIFFSIGN definition are the"
                             " same")
            raise SyntaxError(msg)

        if (r1, r2) in self.samesign_pairs or (r2, r1) in self.samesign_pairs:
            msg = self.error(line_no, "DIFFSIGN definition for reactions "
                             "previously defined as SAMESIGN. Error")
            raise SyntaxError(msg)

        if (r1, r2) in self.diffsign_pairs or (r2, r1) in self.diffsign_pairs:
            msg = "Warning: Duplicate DIFFSIGN definition in line %u" % line_no
            self.messages.append((Verbosity.WARNING, msg))
        else:
            self.diffsign_pairs.add((r1, r2))

    @staticmethod
    def _simplifyLinExpr(lexpr):
        """ simplify the given linear expression by combining coefficients

        This function combines all entries with the same name by adding up the
        corresponding coefficients, preserving the order of the names (first
        occurrence). Terms with coefficient zero are removed.

        Keyword arguments:

        lexpr -- list of (coefficient, name) pairs

        Returns: list of (coefficient, name) pairs with unique name
        """
        firstIndex = {}
        nUniqueNames = 0
        result = []
        for coef, name in lexpr:
            if name not in firstIndex:
                result.append((coef, name))
                firstIndex[name] = nUniqueNames
                nUniqueNames += 1
            else:
                j = firstIndex[name]
                result[j] = ((result[j][0] + coef), name)

        # Remove terms with coefficient zero
        for i in range(len(result)-1, -1, -1):
            if result[i][0] == 0.:
                del result[i]

        return result

    def _simplifyLinIneqs(self):
        """ simplify self.lin_constraints by
            a) combining terms using _simplifyLinExpr
            b) removing terms whose coefficients are 0
            c) converting single-term constraints to LB/UB constraints
        """
        for i in range(len(self.lin_constraints)):
            line_no, line = self.lin_constraints_orig[i]
            lineq = self.lin_constraints[i]

            # First simplify the left-hand side by combining coefficients
            lhs = self._simplifyLinExpr(lineq.lhs)
            nTerms = len(lhs)

            # Make sure that at least one coefficient is not zero
            if nTerms == 0:
                # Mark constraints with all zero coefficients for deletion
                self.lin_constraints[i] = None
                if (lineq.isEq and lineq.rhs != 0.) or (not lineq.isEq and
                                                        lineq.rhs < 0.):
                    msg = ("Warning: Conflicting constraint '%s' in line %u"
                           " - ignored" % (line, line_no))
                    self.messages.append((Verbosity.WARNING, msg))
                else:
                    msg = ("Info: Ignoring tautological constraint '%s' in line"
                           " %u" % (line, line_no))
                    self.messages.append((Verbosity.INFO, msg))
                continue

            # Transform entries with only one term to LB/UB definitions
            if nTerms == 1:
                coef, name = lhs[0]
                lb = self.lb[name] if name in self.lb else -inf
                ub = self.ub[name] if name in self.ub else inf
                tmpValue = old_div(lineq.rhs, coef)

                if lineq.isEq:
                    if name in self.lb:
                        if tmpValue < lb:
                            msg = "'%s' (line %u) violates LB %s %g" % (line,
                                                                        line_no, name, lb)
                            raise SyntaxError(msg)
                        msg = ("Warning: '%s' (line %u) overrides earlier "
                               "LB %s %g" % (line, line_no, name, lb))
                        self.messages.append((Verbosity.WARNING, msg))

                    if name in self.ub:
                        if tmpValue > ub:
                            msg = "'%s' (line %u) violates UB %s %g" % (line,
                                                                        line_no, name, ub)
                            raise SyntaxError(msg)
                        msg = ("Warning: '%s' in line %u overrides earlier "
                               "UB %s %g" % (line, line_no, name, ub))
                        self.messages.append((Verbosity.WARNING, msg))

                    self.lb[name] = self.ub[name] = tmpValue

                elif coef > 0:
                    if name in self.ub:
                        if tmpValue < ub:
                            msg = ("Warning: '%s' (line %u) overrides earlier "
                                   "UB %s %g" % (line, line_no, name, ub))
                            self.messages.append((Verbosity.WARNING, msg))
                            self.ub[name] = tmpValue
                        else:
                            msg = ("Warning: '%s' (line %u) is ignored because "
                                   "it is less strict than earlier UB %s %g" %
                                   (line, line_no, name, ub))
                            self.messages.append((Verbosity.WARNING, msg))
                    else:
                        self.ub[name] = tmpValue

                else:
                    if name in self.lb:
                        if tmpValue < lb:
                            msg = ("Warning: '%s' (line %u) overrides earlier "
                                   "LB %s %g" % (line, line_no, name, lb))
                            self.messages.append((Verbosity.WARNING, msg))
                            self.lb[name] = tmpValue
                        else:
                            msg = ("Warning: '%s' (line %u) is ignored because "
                                   "it is less strict than earlier LB %s %g" %
                                   (line, line_no, name, lb))
                            self.messages.append((Verbosity.WARNING, msg))
                    else:
                        self.lb[name] = tmpValue

                # Mark original constraint for deletion
                self.lin_constraints[i] = None
            else:
                self.lin_constraints[i] = LinearConstraint(lhs, lineq.getSign(),
                                                           lineq.rhs)

        # Actually remove all constraints marked as deleted
        for i in range(len(self.lin_constraints)-1, -1, -1):
            if not self.lin_constraints[i]:
                del self.lin_constraints[i]
                del self.lin_constraints_orig[i]

    def checkBounds(self):
        """ checks whether all pairs of lower and upper bounds (specific or
        unspecific as applicable) are sensible (i.e. lower bound <= upper bound)
        """
        for name in self.lb:
            if name in self.ub:
                if self.lb[name] > self.ub[name]:
                    msg = self.error(0, "Bound mismatch: LB(%s) > UB(%s)" %
                                     (name, name))
                    raise SyntaxError(msg)

    def checkUnique(self, keyword, line_no=0):
        """ check whether occurrence of keyword is unique

        Keyword argument:

        keyword          -- the keyword to check
        msg              -- as in self.error()

        Modified member variables:

        got[keyword]     -- changed from False to True at first occurrence
        """
        if self.got[keyword]:
            msg = self.error(line_no, "Duplicate %s line encountered" %
                             keyword)
            raise SyntaxError(msg)
        self.got[keyword] = True

    def error(self, line_no=0, msg="Syntax error", hint=""):
        """ generate an error message

        Keyword arguments:

        line_no     -- the line number, 0 if not applicable
        msg         -- error message
        hint        -- hint on correct syntax
        """
        s = msg
        if line_no > 0:
            s += " in line %u" % line_no
        if hint != "":
            s += "\nUsage:\n"
            s += hint
        return s
