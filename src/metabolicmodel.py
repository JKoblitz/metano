"""
This module defines the class MetabolicModel, which stores a stoichiometric
metabolic model.


This file is part of metano.
Copyright (C) 2010-2019 Alexander Riemer, Julia Helmecke
Braunschweig University of Technology,
Dept. of Bioinformatics and Biochemistry

Rene Rex contributed to this file.

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

from future import standard_library
standard_library.install_aliases()
from builtins import map
from builtins import zip
from builtins import range
from builtins import object
from metano.defines import Verbosity, ReaFileStruc, typename, makeUnique
from numpy import inf, array, nonzero
from io import StringIO
from metano.reactionparser import ReactionParser
import re
import codecs


class ModelError(Exception):
    """ Exception raised if an inconsistency in the model is detected
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Reactant(object):
    """ Class for a metabolite occurring in a reaction, stored with coefficient

        Member variables:

        name -- the metabolite's name
        coef -- the stoichiometric coefficient of the metabolite in the reaction
                (< 0: consumed, > 0: produced)
    """
    def __init__(self, name="", coef=0.):
        self.name = name
        self.coef = coef

    @staticmethod
    def reactantFromIterable(x):
        """ return a reactant from iterable x: coef <-- x[0], name <-- x[1];
            if x is already a Reactant, return x
        """
        if isinstance(x, Reactant):
            return x
        try:
            return Reactant(*x[1::-1])

        except TypeError:
            raise TypeError("Expected iterable or Reactant. Got '%s'" %
                            typename(x))
        except ValueError:
            raise TypeError("Expected iterable with at least two "
                            "elements")

    def __repr__(self):
        return "Reactant(%r, %r)" % (self.name, self.coef)

    def __iter__(self):
        """ return an iterator over the tuple (coef, name) - for conversion to
            tuple
        """
        return iter((self.coef, self.name))


class Reaction(object):
    """ Represents a single chemical reaction

        Member variables:

        name        -- the reaction's identifier
        direction   -- 0: reversible, 1: left to right, -1: right to left
        metabolites -- list of Reactant objects with name and coefficient
                       (must never be changed directly!); no two elements have
                       the same name (coef < 0: left-hand side,
                                      coef > 0: right-hand side)
        lb          -- lower bound for flux (-inf means unbounded) (property)
        ub          -- upper bound for flux (inf means unbounded)  (property)
        comment     -- comment string, single- or multi-line

        _metabDict  -- dictionary {name : index in metabolites}
        _lb, _ub    -- variables for lower and upper flux bound
    """

    # Mapping from arrow to corresponding direction enum value
    arrowToDirection = {ReaFileStruc.arrowRev: 0, ReaFileStruc.arrowIrr: 1,
                        ReaFileStruc.arrowFlip: -1}
    directionToArrow = dict(list(zip(list(arrowToDirection.values()),
                                list(arrowToDirection.keys()))))

    def __init__(self, name="", educts=None, arrow=None, products=None,
                 comment="", educt_factor=-1., lb=-inf, ub=inf, **kwargs):
        """ create a new reaction

        Keyword arguments:

        name      -- the reaction's identifier
        educts    -- reactants on left-hand side (list of Reactant objects or
                     list of (coef, name) pairs)
        arrow     -- <=>/-->/<--, specifies direction and reversibility
        products  -- reactants on right-hand side (list of Reactant objects or
                     list of (coef, name) pairs))
        educt_factor
                  -- factor for educts coefficients (default -1. assumes that
                     educts have positive signs)
        comment   -- comment string
        lb        -- lower bound for flux (-inf means unbounded)
        ub        -- upper bound for flux (inf means unbounded)

        direction -- reaction direction (alternatively to arrow)
        """
        self.name, self.comment, self._lb, self._ub = name, comment, lb, ub

        if not arrow:
            self.setDirection(kwargs.get("direction", 0))
            # If no direction or arrow is given, reaction is assumed reversible
        else:
            try:
                self.setDirection(self.arrowToDirection[arrow])
                # Arrow takes precedence over a possible "direction" keyword
                # argument
            except KeyError:
                raise ValueError("Illegal reaction arrow '%s'. Must be one of "
                                 "%s" % (arrow, list(self.arrowToDirection.keys())))

        self.metabolites = []
        self._metabDict = {}
        # Convert educts and products to Reactants if given as iterables
        for e in map(Reactant.reactantFromIterable, educts):
            self.addReactant(e, educt_factor)
        for p in map(Reactant.reactantFromIterable, products):
            self.addReactant(p)


    def clear(self):
        """ clear the Reaction object
        """
        self.name = ""
        self.direction = 0
        self.metabolites = []
        self._metabDict = {}
        self.comment = ""
        self._lb = -inf
        self._ub =  inf

    def addReactant(self, reactant, factor=1.):
        """ add the given reactant to the reaction, with its coefficient
            multiplied by the given factor
        """
        # If metabolite is already present, combine reactants (add coefficients)
        if reactant.name in self._metabDict:
            self.metabolites[self._metabDict[reactant.name]].coef += \
                reactant.coef*factor
        else:
            new_index = len(self.metabolites)
            self._metabDict[reactant.name] = new_index
            self.metabolites.append(reactant)
            self.metabolites[new_index].coef *= factor

    def removeReactant(self, name):
        """ remove the reactant given by name from the reaction
        """
        if name not in self._metabDict:
            raise KeyError(name+" not found in reaction "+self.name)

        del self.metabolites[self._metabDict[name]]
        # Update indices in dictionary
        for i in range(self._metabDict[name], len(self.metabolites)):
            self._metabDict[self.metabolites[i].name] = i
        del self._metabDict[name]


    def removeZeroReactants(self, tolerance=0.):
        """ remove all reactants with coefficients that are less than the given
            tolerance from zero
        """
        i = len(self.metabolites)-1
        while i >= 0:
            coef, name = self.metabolites[i]
            if abs(coef) <= tolerance:
                self.removeReactant(name)
            i -= 1


    def setDirection(self, direction):
        """ set reaction direction (either as arrow or as direction enum)
        """
        if direction in self.arrowToDirection:
            self.direction = self.arrowToDirection[direction]
        elif direction in self.directionToArrow:
            self.direction = direction
        else:
            raise ValueError("Illegal direction value or arrow: %r" % direction)

        # Adjust bounds if necessary
        if self.direction == 1:
            self._lb = max(self._lb, 0.)
            self._ub = max(self._ub, 0.)
        elif self.direction == -1:
            self._lb = min(self._lb, 0.)
            self._ub = min(self._ub, 0.)


    def getEducts(self, withZero=False):
        """ return the educts (i.e. reactants with coefficient < 0) as Reactant
            objects; if withZero is True, include reactants with coefficient 0
        """
        if withZero:
            return [x for x in self.metabolites if x.coef <= 0.]
        else:
            return [x for x in self.metabolites if x.coef < 0.]

    def getProducts(self, withZero=False):
        """ return the products (i.e. reactants with coefficient > 0);
            if withZero is True, include reactants with coefficient 0
        """
        if withZero:
            return [x for x in self.metabolites if x.coef >= 0.]
        else:
            return [x for x in self.metabolites if x.coef > 0.]

    def setBounds(self, lb, ub):
        self.lb = lb
        self.ub = ub

    def getBounds(self):
        return self._lb, self._ub

    def getLb(self):
        """ lower bound for the flux through the reaction
        """
        return self._lb

    def setLb(self, val):
        """ set lower bound for flux, adjusted for reaction direction
        """
        if self.direction == 0:
            self._lb = val
        elif self.direction == 1:
            self._lb = max(val, 0.)
        else:
            self._lb = min(val, 0.)

    def getUb(self):
        """ upper bound for the flux through the reaction
        """
        return self._ub

    def setUb(self, val):
        """ set lower bound for flux, adjusted for reaction direction
        """
        if self.direction == 0:
            self._ub = val
        elif self.direction == 1:
            self._ub = max(val, 0.)
        else:
            self._ub = min(val, 0.)

    lb = property(getLb, setLb)
    ub = property(getUb, setUb)


    def writeToFileHandle(self, f, reactionLabels=None, metaboliteLabels=None,
                          width=0, addComment=False, addNewLine=True):
        """ write the reaction equation to the given file handle, replace
            reaction name and reactant names with labels if given

            Keyword arguments:

            f                -- file object (must be open for writing)
            reactionLabels   -- either a single label or dict {name : label}
            metaboliteLabels -- dict {name : label} for metabolites
            width            -- string width for alignment of reaction name
            addComment       -- if True, append '  # <comment string>' to output
            addNewLine       -- if True, end output with newline character
        """
        if reactionLabels:
            try:
                reaname = reactionLabels[self.name]
            except TypeError:
                reaname = reactionLabels
        else:
            reaname = self.name

        s = ""
        if addComment and self.comment and self.comment.find('\n') >= 0:
            multiLine = True
            # Display multi-line comment before reaction
            for line in self.comment.split('\n'):
                s += ReaFileStruc.commentSign + " " + line + '\n'
        else:
            multiLine = False

        # Only show reactants with non-zero coefficients
        educts, products = [], []
        if metaboliteLabels:
            for coef, name in self.metabolites:
                if coef < 0.:
                    educts.append("%g %s" % (-coef,
                                  metaboliteLabels.get(name, name)))
                elif coef > 0.:
                    products.append("%g %s" % (coef,
                                    metaboliteLabels.get(name, name)))
        else:
            for coef, name in self.metabolites:
                if coef < 0.:
                    educts.append("%g %s" % (-coef, name))
                elif coef > 0.:
                    products.append("%g %s" % (coef, name))

        arrow = self.directionToArrow[self.direction]
        s += "%s %s %s %s %s" % (reaname.ljust(width), ReaFileStruc.delimiter,
                                " + ".join(educts), arrow, " + ".join(products))
        if addComment and not multiLine and self.comment:
            s += "  %s %s" % (ReaFileStruc.commentSign, self.comment)
        if addNewLine:
            s += '\n'
        f.write(s)

    def toString(self, includeZeroReactants=False):
        """ generate a string representation of the reaction for output

        If includeZeroReactants is False, reactants with coefficient zero are
        omitted.
        """
        educts, products = [], []
        for coef, name in self.metabolites:
            # Only display reactants with non-zero coefficients
            if coef < 0.:
                educts.append("%g %s" % (-coef, name))
            elif coef > 0.:
                products.append("%g %s" % (coef, name))
            elif includeZeroReactants:
                educts.append("%g %s" % (0., name))
        arrow = self.directionToArrow[self.direction]
        return (self.name + " "+ReaFileStruc.delimiter+" " + " + ".join(educts)
                + " "+arrow+" " + " + ".join(products))


    def __str__(self):
        return self.toString()

    def __repr__(self):
        educts = [(-coef, name) for (coef, name) in self.getEducts(True)]
        products = list(map(tuple, self.getProducts(False)))
        return ("Reaction(%r, %r, %r, %r, %r, %r, %r, %r)" %
                (self.name, educts, self.directionToArrow[self.direction],
                 products, self.comment, -1., self._lb, self._ub))

    def __iter__(self):
        return iter(self.metabolites)

    def __len__(self):
        return len(self.metabolites)


class MetabolicModel(object):
    """ Represents a metabolic model as list of reactions and keeps track of all
        metabolites in the model.

    Member variables:

    reactions       -- a list of Reaction objects
    reactionDict    -- dictionary {name : index in reactions list}
    metabolites     -- a list of metabolite names
    metaboliteDict  -- dictionary {name : index in metabolites list}

    biomass_name        -- name of biomass reaction (None - not set)
    cofactors           -- list of cofactor metabolites
    _reactionLabels     -- dict {reaction name : label} for displaying the model
    _metaboliteLabels   -- dict {metabolite name : label} for displaying the
                           model
    useReactionLabels   -- if True, use _reactionLabels instead of reaction
                           names when displaying the model
    useMetaboliteLabels -- if True, use _metaboliteLabels instead of metabolite
                           names when displaying the model
    showComments        -- if True, display comments when displaying the model
                           (False by default)
    """

    def __init__(self, cofactors=set(), reactionLabels={}, metaboliteLabels={}):
        """ create a new MetabolicModel

        Keyword arguments:

        cofactors        -- a set of metabolites marked as cofactors
        reactionLabels   -- a dictionary of labels to be used in place of the
                            reaction names when displaying the model
        metaboliteLabels -- a dictionary of labels to be used in place of the
                            metabolite names when displaying the model

        Use of labels can be triggered at any time via the Boolean flags
        self.useReactionLabels and self.useMetaboliteLabels.
        """
        self.clearReactions()
        self.cofactors = cofactors
        self.setReactionLabels(reactionLabels)
        self.setMetaboliteLabels(metaboliteLabels)
        self.showComments = False


    def setReactionLabels(self, reactionLabels):
        """ set the reaction labels, which are used in place of the reaction
            names when displaying the model (given as dictionary {name : label})
        """
        self._reactionLabels = reactionLabels
        self.useReactionLabels = bool(reactionLabels)

    def setMetaboliteLabels(self, metaboliteLabels):
        """ set the metabolite labels, which are used in place of the metabolite
            names when displaying the model (given as dictionary {name : label})
        """
        self._metaboliteLabels = metaboliteLabels
        self.useMetaboliteLabels = bool(metaboliteLabels)

    def makeLabelsUnique(self, doReactionLabels=True, doMetaboliteLabels=True):
        """ generate unique labels for reactions, metabolites, or both
            - identical labels are distinguished by appending numbers
        """
        if doReactionLabels:
            reanames = self.getReactionNames()
            if len(self.reactions) != len(set(reanames)):
                raise ModelError("Reaction names are not unique")
            self._reactionLabels = dict(list(zip(reanames,
                makeUnique(self._getTranslatedReactionNames()))))

        if doMetaboliteLabels:
            # Metabolite names are unique by construction

            self._metaboliteLabels = dict(list(zip(self.getMetaboliteNames(),
                makeUnique(self._getTranslatedMetaboliteNames()))))


    def getBounds(self):
        """ return the bounds of all reaction fluxes as a pair of lists, indexed
            like the list of reactions (and thus like the columns of the
            stoichiometric matrix); -inf/inf means that a reaction is unbounded
            to the left/right side
        """
        lbVec, ubVec = [], []
        for r in self.reactions:
            lbVec.append(r.lb)
            ubVec.append(r.ub)
        return lbVec, ubVec

    def setBounds(self, lbVec=None, ubVec=None):
        """ set lower and upper bounds of all reaction fluxes - both lists are
            indexed like the list of reactions (and thus like the columns of the
            stoichiometric matrix); -inf/inf means that a reaction is unbounded
            to the left/right side

        Keyword arguments:

        lbVec -- vector of lower bounds (may be empty)
        ubVec -- vector of upper bounds (may be empty)
        """
        nReactions = len(self.reactions)
        if lbVec is not None:
            if len(lbVec) != nReactions:
                raise ValueError("lbVec has wrong length (%u); expected %u" %
                                 (len(lbVec), nReactions))
            for i in range(nReactions):
                self.reactions[i].lb = lbVec[i]
        if ubVec is not None:
            if len(ubVec) != nReactions:
                raise ValueError("ubVec has wrong length (%u); expected %u" %
                                 (len(ubVec), nReactions))
            for i in range(nReactions):
                self.reactions[i].ub = ubVec[i]

    def getFiniteBounds(self):
        """ return the finite bounds of all reaction fluxes as a pair of dicts
            {reaction name : value}
        """
        lbDict, ubDict = {}, {}
        for r in self.reactions:
            lb, ub = r.lb, r.ub
            if lb > -inf:
                lbDict[r.name] = lb
            if ub < inf:
                ubDict[r.name] = ub
        return lbDict, ubDict

    def setFiniteBounds(self, lbDict={}, ubDict={}, strict=True, messages=[]):
        """ set lower and upper bounds of fluxes through specific reactions

        Keyword arguments:

        lbDict   -- dict {reaction name : value} for lower bounds
        ubDict   -- dict {reaction name : value} for upper bounds
        strict   -- if True, raise ValueError if both bounds are outside the
                    range allowed by the reaction's direction
        messages -- list in which to store warnings and info messages (as
                    (level, msg) pairs)
        """
        for reaname in ubDict:
            try:
                index = self.reactionDict[reaname]
            except KeyError:
                msg = "Warning: UB set for non-existing reaction %s" % reaname
                messages.append((Verbosity.WARNING, msg))
                continue

            if self.reactions[index].direction == 1 and ubDict[reaname] < 0.:
                msg = ("Error: Upper bound for irreversible reaction %s is "
                       "negative" % reaname)
                if strict:
                    raise ValueError(msg)
                else:
                    messages.append((Verbosity.ERROR, msg))
            self.reactions[index].ub = ubDict[reaname]

        for reaname in lbDict:
            try:
                index = self.reactionDict[reaname]
            except KeyError:
                msg = "Warning: LB set for non-existing reaction %s" % reaname
                messages.append((Verbosity.WARNING, msg))
                continue

            if self.reactions[index].direction == -1 and lbDict[reaname] > 0.:
                msg = ("Error: Lower bound for right-to-left irreversible "
                       "reaction %s is negative" % reaname)
                if strict:
                    raise ValueError(msg)
                else:
                    messages.append((Verbosity.ERROR, msg))
            self.reactions[index].lb = lbDict[reaname]


    def getStoichiometricMatrix(self):
        """ compute the stoichiometric matrix from the list of reactions

        Rows are indexed like the metabolite list, columns are indexed like the
        reaction list.
        """
        nReactions = len(self.reactions)
        matrix = [[0.]*nReactions for _ in self.metaboliteDict]
        for j in range(nReactions):
            for coef, name in self.reactions[j]:
                matrix[self.metaboliteDict[name]][j] = coef

        return matrix


    def clearReactions(self):
        """ clear the network, i.e. remove all reactions and metabolites
        """
        self.reactions = []
        self.reactionDict = {}
        self.metabolites = []
        self.metaboliteDict = {}
        self.biomass_name = None


    def addReactions(self, reactions):
        """ append the given list of reactions to the model

            Reactions are given in one of the following formats:
            a) Reaction objects
            b) tuples (name, educts, arrow, products[, comment]) with educts and
               products either as Reactant objects or (coef, name) tuples
            c) tuples (name, educts, arrow, products, comment, lb, ub)
            In cases b) and c), educts must have positive coefficients.
        """
        first_new = len(self.reactions)
        reactionSet = set(self.reactionDict.keys())
        for r in reactions:
            if isinstance(r, Reaction):
                if r.name in reactionSet:
                    raise ModelError("Reaction '%s' already in model" % r.name)
                reactionSet.add(r.name)
                self.reactions.append(r)
            else:
                try:
                    numElements = len(r)
                except TypeError:
                    raise TypeError("Expected tuple or Reaction, got '%s'" %
                                    typename(r))
                if numElements < 4:
                    raise TypeError("Reaction tuple has wrong length - must be "
                                    "at least (name, educts, arrow, products)")
                reaname = r[0]
                if reaname in reactionSet:
                    raise ModelError("Reaction '%s' already in model" % reaname)
                reactionSet.add(reaname)
                if numElements < 6:
                    self.reactions.append(Reaction(*r[:5]))
                else:
                    self.reactions.append(Reaction(*(r[:5]+(-1.,)+r[5:7])))

            # If any reaction has name 'biomass' (case-independent), mark it as
            # biomass reaction if biomass_name is not set otherwise
            if (self.reactions[-1].name.lower() == "biomass" and
                not self.biomass_name):
                self.biomass_name = self.reactions[-1].name

        # Update reactionDict, metabolites, and metaboliteDict
        nMetabolites = len(self.metabolites)
        for i in range(first_new, len(self.reactions)):
            reaction = self.reactions[i]
            self.reactionDict[reaction.name] = i

            for reactant in reaction:
                met = reactant.name
                if met not in self.metaboliteDict:
                    self.metaboliteDict[met] = nMetabolites
                    self.metabolites.append(met)
                    nMetabolites += 1


    def addReaction(self, reaction):
        """ add a single reaction to the model using addReactions - see
            definition of addReactions for a description of the input format
        """
        self.addReactions([reaction])


    def addReactionsFromFile(self, filename, rparser=None):
        """ parse the file given by filename; if no ReactionParser is supplied,
            one is created for the occasion
        """
        if not rparser:
            rparser = ReactionParser()
        self.addReactions(rparser.parse(filename))


    def addReactionsFromFileHandle(self, f, rparser=None):
        """ parse the file given by file handle f
        """
        if not rparser:
            rparser = ReactionParser()
        self.addReactions(rparser.parseByHandle(f))

    def addReactionsFromString(self, s, rparser=None):
        """ parse the reaction file contents given in string s
        """
        if not rparser:
            rparser = ReactionParser()
        reafile = StringIO(s+'\n')
        try:
            reactions = rparser.parseByHandle(reafile)
        finally:
            reafile.close()
        self.addReactions(reactions)

    def addReactionsFromList(self, reactionList, rparser=None):
        """ parse the given list of reaction strings
        """
        self.addReactionsFromString('\n'.join(reactionList), rparser)


    def removeReaction(self, r, doUpdate=True):
        """ remove the reaction given by name, by index, or as Reaction object
            from the list; metabolites and metaboliteDict are only updated if
            doUpdate is True
        """
        if isinstance(r, Reaction):
            name = r.name
            if name not in self.reactionDict:
                raise KeyError(name+" not found in model")
            index = self.reactionDict[name]
        elif isinstance(r, str):
            name = r
            if name not in self.reactionDict:
                raise KeyError(name+" not found in model")
            index = self.reactionDict[name]
        else:
            index = r
            try:
                name = self.reactions[index].name
            except TypeError:
                raise TypeError("Expected Reaction, str, or int, got '%s'" %
                                typename(r))

        del self.reactions[index]
        del self.reactionDict[name]
        # Update reactionDict
        for i in range(index, len(self.reactions)):
            self.reactionDict[self.reactions[i].name] = i

        if doUpdate:
            self._updateMetaboliteList()


    def _updateMetaboliteList(self):
        """ build metabolites list and metaboliteDict from reactions list

        This function is called internally by removeReaction(). It preserves the
        order of the metabolites list. However, a call to removeReaction()
        followed by addReactions() will most likely change the order.
        Note: This function only accesses the reactions list, not reactionDict.
        """
        # Set of all metabolites in all reactions (for checking)
        metabSet = set()
        # List of all metabolites in all reactions (for order)
        metabolites = []
        for reaction in self.reactions:
            for _, met in reaction:  # get only names, ignore coefficients
                if met not in metabSet:
                    metabSet.add(met)
                    metabolites.append(met)

        # Remove all metabolites from model that are not present in metabSet
        nMetabolites = len(self.metabolites)
        first_deleted = nMetabolites
        for i in range(nMetabolites-1, -1, -1):
            met = self.metabolites[i]
            if met not in metabSet:
                del self.metabolites[i]
                del self.metaboliteDict[met]
                first_deleted = i

        # Update self.metaboliteDict
        nMetabolites = len(self.metabolites)
        for i in range(first_deleted, nMetabolites):
            self.metaboliteDict[self.metabolites[i]] = i
        # As this function is called only after removing reactions, we don't
        # need to check for metabolites in list 'metabolites' that are not
        # present in self.metabolites.


    def removeReactionsByNameList(self, nameList, doUpdate=True):
        """ remove all reactions specified by name in nameList;
            metabolites and metaboliteDict are only updated if doUpdate is True

        Names that don't correspond to any reaction in the model are ignored.
        """
        for name in nameList:
            if name in self.reactionDict:
                self.removeReaction(name, False)
        if doUpdate:
            self._updateMetaboliteList()

    def removeReactionsByList(self, reactionList, doUpdate=True):
        """ remove all reactions specified by the strings (reaction equations)
            in reactionList;
            metabolites and metaboliteDict are only updated if doUpdate is True

        Reactions are identified by name only, equations are not checked.
        """
        rparser = ReactionParser()
        reafile = StringIO('\n'.join(reactionList)+'\n')
        try:
            for reaction_entry in rparser.parseByHandle(reafile):
                self.removeReaction(reaction_entry[0], False)
        finally:
            reafile.close()
        if doUpdate:
            self._updateMetaboliteList()

    def removeReactionsByModel(self, model, doUpdate=True):
        """ remove all reactions present in both model and self from self;
            metabolites and metaboliteDict are only updated if doUpdate is True

        Reactions are identified by name only, Reactants are not checked.
        """
        for reaction in model:
            if reaction.name in self.reactionDict:
                self.removeReaction(reaction, False)
        if doUpdate:
            self._updateMetaboliteList()


    def getSubmodelByList(self, nameList):
        """ construct a submodel from the reactions identified by nameList

        Names that don't correspond to any reaction in the model are ignored.
        """
        m = MetabolicModel(self.cofactors, self._reactionLabels,
                           self._metaboliteLabels)
        # Collect all reactions whose names are in nameList
        rlist = []
        for name in nameList:
            if name in self.reactionDict:
                rlist.append(self.reactions[self.reactionDict[name]])
        m.addReactions(rlist)
        return m


    def getSubModelByExcludeList(self, nameList):
        """ construct a submodel from all reactions in the model not in nameList

        Names that don't correspond to any reaction in the model are ignored.
        """
        m = MetabolicModel(self.cofactors, self._reactionLabels,
                           self._metaboliteLabels)
        # Collect all reactions whose names are not in nameList
        nameSet = set(nameList)
        rlist = []
        for rea in self.reactions:
            if rea.name not in nameSet:
                rlist.append(rea)
        m.addReactions(rlist)
        return m


    def getSubmodelByRegex(self, pattern):
        """ construct submodel from the reactions identified by the name pattern

        Pattern can be either a single regular expression (given as string) or a
        list of regular expressions.
        """
        if isinstance(pattern, str):
            pattern = [pattern]

        regex = "|".join("^.*"+p+".*$" for p in pattern)
        matchedReactions = re.findall(re.compile(regex, re.MULTILINE),
                                        self.writeToString(False))

        m = MetabolicModel(self.cofactors, self._reactionLabels,
                           self._metaboliteLabels)
        m.addReactionsFromList(matchedReactions)
        # Add bounds and comments lost due to exporting and reading reactions
        # from string
        if self.useReactionLabels:
            # Create reverse mapping of labels to reaction names
            revMap = dict(list(zip(list(self._reactionLabels.values()),
                              list(self._reactionLabels.keys()))))
        else:
            rNames = self.getReactionNames()
            revMap = dict(list(zip(rNames, rNames)))
        for rea in m.reactions:
            original = self.reactions[self.reactionDict[revMap[rea.name]]]
            rea.setBounds(*original.getBounds())
            rea.comment = original.comment

        return m


    def getReactionNames(self):
        """ return the names (not labels) of the reactions in the model (in
            correct order)
        """
        return [r.name for r in self.reactions]

    def getMetaboliteNames(self):
        """ return the names (not labels) of the metabolites in the model (in
            correct order)
        """
        return self.metabolites


    def compare(self, otherModel):
        diffModel = MetabolicModel()
        # Compute difference based on line-by-line comparison of exported models
        self_export = self.writeToString(False).split('\n')
        other_export = otherModel.writeToString(False).split('\n')
        diffModel.addReactionsFromList(list(
            set(self_export).difference(other_export)))
        return diffModel


    def _getTranslatedReactionNames(self):
        if self.useReactionLabels:
            return [self._reactionLabels.get(reaction.name,
                                             reaction.name).replace(" ", "_")
                    for reaction in self.reactions]
        else:
            return self.getReactionNames()

    def _getTranslatedMetaboliteNames(self):
        if self.useMetaboliteLabels:
            return [self._metaboliteLabels.get(met, met).replace(" ", "_")
                    for met in self.metabolites]
        else:
            return self.getMetaboliteNames()

    def writeToFileHandle(self, outputFile, withComments=True,
                          useReactionLabels=None, useMetaboliteLabels=None,
                          sort=False):
        """ write the model (as reaction file) to the given file handle

        Keyword arguments:

        outputFile          -- file object (must be open for writing)
        withComments        -- if True, include reaction comments
        useReactionLabels   -- if not None, self.~ is set to the given value
        useMetaboliteLabels -- if not None, self.~ is set to the given value
        sort                -- if True, the reactions will be sorted by reaction
                               name
        """
        if not self.reactions:
            return  # Nothing to do

        if useReactionLabels is not None:
            self.useReactionLabels = useReactionLabels
        if useMetaboliteLabels is not None:
            self.useMetaboliteLabels = useMetaboliteLabels

        reactionLabels = self._getTranslatedReactionNames()
        maxlen = len(max(reactionLabels, key=len))
        rlabels = dict(list(zip(self.getReactionNames(), reactionLabels)))
        mlabels = dict(list(zip(self.getMetaboliteNames(),
                           self._getTranslatedMetaboliteNames())))

        last_index = len(self.reactions)-1
        count = 0

        if sort:

            for reactionName in sorted(self.getReactionNames()):
                reaction = self.reactions[self.reactionDict[reactionName]]
                reaction.writeToFileHandle(outputFile, rlabels, mlabels, maxlen,
                                           withComments, count < last_index)
                count += 1
        else:

            for reaction in self.reactions:
                reaction.writeToFileHandle(outputFile, rlabels, mlabels, maxlen,
                                           withComments, count < last_index)
                count += 1

    def writeToFile(self, outputFilename, withComments=True,
                    useReactionLabels=None, useMetaboliteLabels=None,
                    sort=False):
        """ write the model (as reaction file) to the given file using
            self.writeToFileHandle()
        """
        with codecs.open(outputFilename, 'w', encoding='utf-8') as outputFile:
            self.writeToFileHandle(outputFile, withComments, useReactionLabels,
                                   useMetaboliteLabels, sort)
            outputFile.write('\n')

    def writeToString(self, withComments=None, useReactionLabels=None,
                      useMetaboliteLabels=None):
        """ return a string representation of the model (as reaction file);

        Keyword arguments:

        withComments        -- if True, include reaction comments;
                               if None, use self.showComments for decision
        useReactionLabels   -- if not None, self.~ is set to the given value
        useMetaboliteLabels -- if not None, self.~ is set to the given value
        """
        outputFile = StringIO()
        if withComments is None:
            with_comments = self.showComments
        else:
            with_comments = withComments

        self.writeToFileHandle(outputFile, with_comments, useReactionLabels,
                               useMetaboliteLabels)

        output = outputFile.getvalue()
        outputFile.close()
        return output

    def __str__(self):
        return self.writeToString()


    def __iter__(self):
        return iter(self.reactions)

    def __len__(self):
        return len(self.reactions)


    def writeBoundsToParamFileHandle(self, f, useReactionLabels=None):
        """ write all bounds as "LB|UB <reaction name> <value>" lines to the
            given file handle (expected to be open for writing)

        Keyword arguments:

        f                   -- file object (must be open for writing)
        useReactionLabels   -- if not None, self.~ is set to the given value
        """
        if not self.reactions:
            return  # Nothing to do

        if useReactionLabels is not None:
            self.useReactionLabels = useReactionLabels

        reactionLabels = self._getTranslatedReactionNames()
        rlabels = dict(list(zip(self.getReactionNames(), reactionLabels)))

        for rea in self.reactions:
            lb, ub = rea.getBounds()
            if lb > -inf:
                f.write("LB %s %r\n" % (rlabels[rea.name], lb))
            if ub < inf:
                f.write("UB %s %r\n" % (rlabels[rea.name], ub))

    # ------ Exporters for different formats -----------------------------------

    def _formatReactantAnnetCsv(self, coef, name, is_ext, epsilon=1E-8):
        """ format a single reactant given with coefficient for output in anNET
            format (called internally by toAnnetCsvByHandle())

        Keyword arguments:
        coef    -- coefficient
        name    -- reactant's name
        is_ext  -- if True, mark as external, if False mark as cytosolic,
                   otherwise don't mark
        epsilon -- minimum float to be considered different from zero
        """
        # If coefficient is 1, omit it
        if abs(coef-1.) >= epsilon:
            s = "(%r) %s" % (coef, name)
        else:
            s = name
        try:
            c = {True : 'e', False : 'c'}[is_ext]
            s += "[%c]" % c
        except KeyError:
            pass
        return s

    def toAnnetCsv(self, filename, epsilon=1E-8):
        """ export the model as CSV file in anNET format to the file given by
            filename
        """
        with open(filename, 'w') as f:
            self.toAnnetCsvByHandle(f, epsilon)

    def toAnnetCsvByHandle(self, f, epsilon=1E-8):
        """ export the model as CSV file in anNET format to the given file

        Keyword arguments:
        f         -- file object (must be open for writing)
        epsilon   -- threshold below which a float is to be considered zero;
                     must be positive!
        """
        nMetabolites = len(self.metabolites)
        if nMetabolites == 0:
            return      # nothing to do

        # Replace parentheses with underscores in metabolite names
        metabReplace = dict(list(zip(self.metabolites,
                                [met.replace('(', '_').replace(')', '_')
                                 for met in self.metabolites])))

        # Write E. coli compartment data (must be adapted manually)
        f.write(";ID;pH;IS;Potential mV;Volume;\n"
                "compartment;e;3.5;0.15;0;0;extracellular\n"
                "compartment;c;6.5;0.15;0;1;cytosol\n\n;Model;;;;;\n\n"
                ";Abbreviation;reactions;;;;\n")
        # Write reactions
        maxlen = len(max(self.reactions, key=lambda x : len(x.name)))
        for rea in self.reactions:
            arrow = Reaction.directionToArrow[rea.direction]
            # Get educts and products as (coef, name, is_external) tuples
            is_external = lambda x : x[-3:].lower() == "_ex"
            educts = [(-x.coef, metabReplace[x.name], is_external(x.name))
                      for x in rea.getEducts()]
            products = [(x.coef, metabReplace[x.name], is_external(x.name))
                        for x in rea.getProducts()]

            if rea.name == self.biomass_name:
                # Save biomass for last
                biomass_data = (educts, products, arrow)
                continue

            line = "reaction;"+rea.name.ljust(maxlen)+";"
            any_ext = False
            for (_, _, is_ext) in educts+products:
                if is_ext:
                    any_ext = True
                    break
            if not any_ext:
                educts = [(coef, name, None) for (coef, name, _) in educts]
                products = [(coef, name, None) for (coef, name, _) in products]
                line += "[c]"  # cytosolic reaction
            line += (" + ".join([self._formatReactantAnnetCsv(*e+(epsilon,))
                                for e in educts]) + " %s " % arrow +
                     " + ".join([self._formatReactantAnnetCsv(*p+(epsilon,))
                                for p in products]))
            f.write(line+'\n')

        educts, products, arrow = biomass_data
        line = ("\n;Biomass Reaction;\nreaction;%s;" % self.biomass_name +
                " + ".join([self._formatReactantAnnetCsv(e[0], e[1], None,
                            epsilon) for e in educts]) + " %s " % arrow +
                " + ".join([self._formatReactantAnnetCsv(p[0], p[1], None,
                            epsilon) for p in products]))
        f.write(line+"\n\n\nThermo names;;\n\n")

    # ------ Algorithms --------------------------------------------------------

    def getBoundaryFluxes(self):
        """ identify boundary fluxes

        Returns:
        boundaryFluxes -- list of reactions which either have no products or
                          no educts
        """

        return [reaction.name for reaction in self
                if (not reaction.getEducts()) or (not reaction.getProducts())]

    def getInternalFluxes(self):
        """ identify internal reactions

        Returns list of reactions which have at least one educt and at least one
        product
        """

        return [reaction.name for reaction in self
                if reaction.getEducts() and reaction.getProducts()]


    def canBeZero(self, reactions, strict=False):
        """ for each given reaction, return whether the flux can be zero

        Keyword arguments:

        reactions -- list of reaction names
        strict    -- switch determining the behavior for reactions not present
                     in the model (True: raise KeyError, False: mark positive)

        Returns list of Boolean values indexed like reactions list
        """
        result = []
        for reaction in reactions:
            try:
                rea = self.reactions[self.reactionDict[reaction]]
                result.append(rea.lb <= 0. and rea.ub >= 0.)
            except KeyError:
                if strict:
                    raise
                else:
                    result.append(True)
        return result


    def findDuplicateReactions(self, verbose=False, epsilon=1E-8):
        """ identify duplicate reactions

        Two or more reactions in a model may be identical although they have
        different names or the reactants have a different order. Regarding the
        FBA, these reactions are equivalent and can be removed to improve the
        performance.

        Keyword arguments:

        verbose   -- if True, print debug information to the console
        epsilon   -- threshold below which a float is to be considered zero;
                     must be positive!

        Returns:
        duplicateReactions

        duplicateReactions -- set of reaction names which can be removed without
                              altering the capabilities of the model

        example/doctest:

        >>> m = MetabolicModel()
        >>> m.addReactionsFromString('''
        ... Test1 : A <=> B
        ... Test2 : A --> B
        ... Test3 : A <-- B
        ... Test4 : A <=> B
        ... Test5 : B <=> A
        ... Test6 : B <-- A''')
        >>> duplicateReactions = m.findDuplicateReactions()
        >>> m.removeReactionsByNameList(duplicateReactions)
        >>> print m
        Test1 : 1 A <=> 1 B
        Test2 : 1 A --> 1 B
        Test3 : 1 A <-- 1 B

        """

        matrix = self.getStoichiometricMatrix()
        reactionNames = array(self.getReactionNames())

        if len(matrix) <= 0:
            return []

        # add a row with reaction directions
        matrix.append([reaction.direction for reaction in self])

        # list of lists to numpy array
        matrix = array(matrix)
        _, nReactions = matrix.shape

        # make sure two identical reactions have the reactants on the same side
        for i in range(nReactions):
            nonZeros = nonzero( matrix[:,i])[0] # indices of the reactants
            # swap the reactants and the direction if necessary
            if nonZeros.any() and matrix[nonZeros[0],i] > 0:
                matrix[:,i] = matrix[:,i]*-1

        if verbose:
            print("stoichiometric matrix and reaction names:")
            print(array(matrix), "\n", reactionNames)

        # transpose matrix and create a list of lists
        matrix = [[matrix[j][i] for j in range(len(matrix))]
                                for i in range(len(matrix[0]))]

        # sort matrix and reaction names
        reactionNames = sorted(reactionNames,
                                key=lambda x: matrix[self.reactionDict[x]])
        matrix = sorted(matrix)

        if verbose:
            print("transposed and sorted stoichiometric matrix, sorted "\
                    "reaction names:")
            print(array(matrix), "\n", reactionNames)

        matrix = array(matrix)

        # search for duplicates in the sorted matrix
        duplicateReactions = [reactionNames[i] for i in range(1, nReactions)
                                if sum(abs(matrix[i-1,:]-matrix[i,:]))<=epsilon]

        if verbose:
            print("duplicate reactions:")
            print(duplicateReactions)

        return duplicateReactions


    def findDeadEnds(self, recursive=False, unbound=set(), verbose=False,
                     epsilon=1E-8):
        """ identify dead-end metabolites

        A dead end is a metabolite that cannot be produced and consumed in
        separate reactions. Reactions that contain dead-end metabolites always
        get assigned a flux of zero in flux-balance analysis.

        If recursive is True, this function recursively searches for dead ends
        and nonfunctional reactions (which in turn can induce further dead
        ends). Else, it only returns the list of dead-end metabolites.

        Keyword arguments:

        recursive -- if True, work recursively (see explanation above)
        unbound   -- set of reactions to be treated as unbounded (useful for
                     transporters, of which usually only a few are allowed at a
                     time)
        verbose   -- if True, print debug information to the console
        epsilon   -- threshold below which a float is to be considered zero;
                     must be positive!

        Returns:
        deadEnds, deadReactions

        deadEnds      -- list of dead-end metabolites (names)
        deadReactions -- list of reactions containing dead ends
                         (None if recursive is False)
        """
        nMetabolites = len(self.metabolites)
        nReactions = len(self.reactions)

        # First compile lists of sets of producing and consuming reactions for
        # every metabolite and sets of metabolites by reaction in which they
        # participate
        # producingReactions and consumingReactions take bounds into account
        producingReactions = [set() for _ in range(nMetabolites)]
        consumingReactions = [set() for _ in range(nMetabolites)]
        metabolitesByReaction = [set() for _ in range(nReactions)]

        allDeadReactions = set()

        for j in range(nReactions):
            rea = self.reactions[j]
            isUnbounded = rea.name in unbound
            canBePos = rea.ub > epsilon
            canBeNeg = rea.lb < -epsilon
            if not isUnbounded and not canBePos and not canBeNeg:
                # Mark reaction as dead if its flux is constrained to zero
                allDeadReactions.add(j)
                continue

            for coef, name in rea:
                i = self.metaboliteDict[name]

                if abs(coef) > epsilon:
                    metabolitesByReaction[j].add(i)
                    if coef > 0.:
                        # Coefficient is positive, i.e. if the flux is
                        # positive, the metabolite is produced
                        if canBePos or isUnbounded:
                            producingReactions[i].add(j)
                        if canBeNeg or isUnbounded:
                            consumingReactions[i].add(j)
                    else:
                        # Coefficient is negative, i.e. if the flux is
                        # positive, the metabolite is consumed
                        if canBePos or isUnbounded:
                            consumingReactions[i].add(j)
                        if canBeNeg or isUnbounded:
                            producingReactions[i].add(j)

        if verbose:
            nDeadReactions = len(allDeadReactions)
            if nDeadReactions:
                print ("%u reactions have a flux that is restricted to zero:" %
                       nDeadReactions)
                print("\n".join(sorted(self.reactions[i].name for i in
                                       allDeadReactions)))
            else:
                print("There are no reactions with a flux restricted to zero.")
            print()
            print("Producing and consuming reactions for all metabolites:")
            for i in range(nMetabolites):
                print ("%s (%u, %u)" % (self.metabolites[i],
                       len(producingReactions[i]), len(consumingReactions[i])))
            print()

        # Total reactions is the combined set of consuming and producing
        # reactions
        reactionsByMetabolite = [consumingReactions[i] |
            producingReactions[i] for i in range(nMetabolites)]

        isDeadEnd = [False]*nMetabolites

        deadReactions = 1   # dummy value (to enter loop)
        while deadReactions:
            deadReactions = set()

            # 1. Identify dead ends, check for dead reactions

            for i in range(nMetabolites):
                if isDeadEnd[i]:
                    continue  # Don't test known dead ends again

                if (not producingReactions[i] or not consumingReactions[i] or
                    len(reactionsByMetabolite[i]) <= 1):

                    # Metabolite i is a dead end
                    isDeadEnd[i] = True

                    # All reactions that involve this metabolite are also "dead"
                    if verbose:
                        newlyDeadReactions = set()
                        for j in reactionsByMetabolite[i]:
                            deadReactions.add(j)
                            newlyDeadReactions.add(self.reactions[j].name)

                        print("%s prod: %s consu: %s is a dead end (reactions: "
                              "%s)" % (self.metabolites[i],
                                       [self.reactions[j].name
                                        for j in producingReactions[i]],
                                       [self.reactions[j].name for j in
                                        consumingReactions[i]],
                                       ", ".join(newlyDeadReactions)))
                    else:
                        for j in reactionsByMetabolite[i]:
                            deadReactions.add(j)

            if not recursive:
                deadReactions = None
                break

            # 2. Update producingReactions & consumingReactions accordingly

            for j in deadReactions:
                if verbose:
                    print(("Fixing reaction %s (metabolites:" %
                           self.reactions[j].name), end=' ')

                # Remove dead reaction from sets of producing, consuming, and
                # total reactions of every involved metabolite
                for i in metabolitesByReaction[j]:
                    if not isDeadEnd[i]:
                        tmpSet = set((j,))
                        producingReactions[i] -= tmpSet
                        consumingReactions[i] -= tmpSet
                        reactionsByMetabolite[i] -= tmpSet

                allDeadReactions.add(j)

                if verbose:
                    print("%s)" % ", ".join([self.metabolites[i] for i in
                                metabolitesByReaction[j] if not isDeadEnd[i]]))

        # Replace indices with names
        deadEnds = [self.metabolites[i] for i in range(nMetabolites)
                    if isDeadEnd[i]]
        deadReactions = [self.reactions[j].name for j in allDeadReactions]

        if verbose:
            print()
            if recursive:
                print ("Producing and consuming reactions after removal of all "
                       "dead ends:")
                for i in range(nMetabolites):
                    print ("%s (%u, %u)" % (self.metabolites[i],
                        len(producingReactions[i]), len(consumingReactions[i])))
                print()

        return deadEnds, deadReactions


    def getActiveReactions(self, solution, threshold=1E-6, messages=[],
                           coefCutoff=1E-12):
        """ identify the reactions that have an absolute flux above the given
            threshold and all metabolites involved in these reactions

        Keyword arguments:

        solution   -- vector of flux values given as either
                        a) list of values with same order as reactions list
                        b) dictionary {reaction name : value}
        threshold  -- flux threshold above which a reaction is considered
                      "active"
        messages   -- list in which to store warnings and info messages (as
                      (level, msg) pairs)
        coefCutoff -- cutoff below which a stoichiometric coefficient is
                      considered zero

        If solution is given as a dict, all reactions which don't have an entry
        in solution, are considered "inactive". If solution is given as another
        iterable, its length must be the same as the number of reactions in the
        model.

        Returns:
        activeReactions, activeMetabolites

        activeReactions   -- list of active reactions
        activeMetabolites -- dict {name : #occurrences} listing the number of
                             active reactions in which each "active" metabolite
                             occurs
        """
        if isinstance(solution, dict):
            flux = solution
        else:
            if len(solution) != len(self.reactions):
                raise ValueError("Length of solution vector and number of "
                                 "reactions in the model disagree")
            flux = dict(list(zip(self.getReactionNames(), solution)))
        activeMetabolites = {}
        activeReactions = []
        for reaname in flux:
            if reaname not in self.reactionDict:
                msg = ("Warning: Reaction '%s' from solution not found in "
                       "model." % reaname)
                messages.append((Verbosity.WARNING, msg))
                continue

            if abs(flux[reaname]) > threshold:
                activeReactions.append(reaname)
                rea_index = self.reactionDict[reaname]

                for coef, name in self.reactions[rea_index]:
                    if abs(coef) > coefCutoff:
                        if name in activeMetabolites:
                            activeMetabolites[name] += 1
                        else:
                            activeMetabolites[name] = 1

        return activeReactions, activeMetabolites
