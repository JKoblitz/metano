"""
This module defines the class MetabolicFlux, which represents a vector of fluxes
through the reactions of a MetabolicModel.


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
from __future__ import division

from builtins import zip
from builtins import map
from builtins import object
from past.utils import old_div
import csv
from numpy import inf
from metano.defines import typename, addToDictPlus, padNumber
from metano.metabolicmodel import MetabolicModel


class MetabolicFlux(object):
    """ Represents a list of fluxes through a metabolic network. The list is
        stored internally as a dict {reaction name : flux value}.

    Member variables:
    fluxDict   -- dict {name : flux value}
    boundsDict -- dict {name : (lb, ub)} of lower and upper bounds for each flux
    """

    def __init__(self, reactions=None, fluxVec=None, bounds=None):
        """ initialize a metabolic flux instance from one of the following:
            a) (nothing)
            b) dictionary {reaction name : flux value}
            c) MetabolicModel, list of values
            d) list of names, list of values

        Keyword arguments:
        reactions -- dictionary, MetabolicModel, or list of names (see above)
        fluxVec   -- list of flux values (also see above)
        bounds    -- dictionary {reaction name : (lower bound, upper bound)}
        """
        self.clear()
        if not reactions:
            return
        if isinstance(reactions, dict):
            self.fluxDict = reactions
        else:
            try:
                if isinstance(reactions, MetabolicModel):
                    self.setByModel(reactions, fluxVec)
                else:
                    self.setByList(reactions, fluxVec)
            except TypeError:
                raise TypeError("Expected dict or two lists or MetabolicModel "
                                "and list; got (%s, %s)" % (typename(reactions),
                                                            typename(fluxVec)))
        if bounds:
            self.boundsDict = bounds

    def __len__(self):
        return len(self.fluxDict)

    def __getitem__(self, key):
        return self.fluxDict[key]

    def __setitem__(self, key, value):
        self.fluxDict[key] = value

    def __delitem__(self, key):
        del self.fluxDict[key]

    def __iter__(self):
        return iter(self.fluxDict)

    def __repr__(self):
        if not self.fluxDict and not self.boundsDict:
            return "MetabolicFlux()"
        return "MetabolicFlux(%r, bounds=%r)" % (self.fluxDict, self.boundsDict)

    def __str__(self):
        if not self.fluxDict:
            return "MetabolicFlux()"  # Representation of empty object
        s = ""
        # Get length of longest reaction name (for formatting)
        maxlen = len(max(self.fluxDict, key=len))
        nFluxes = len(self.fluxDict)
        count = 0
        for name in sorted(self.fluxDict):
            try:
                lb, ub = self.boundsDict[name]
            except KeyError:
                lb, ub = -inf, inf
            s += "%s : %g (%g, %g)" % (name.ljust(maxlen), self.fluxDict[name],
                                       lb, ub)
            count += 1
            if count < nFluxes:
                s += "\n"
        return s

    def clear(self):
        self.fluxDict = {}
        self.boundsDict = {}

    def clearBounds(self):
        self.boundsDict = {}

    def setBounds(self, name, lb, ub):
        self.boundsDict[name] = lb, ub

    def setToZero(self):
        """ sets all fluxes to zero
        """
        for name in self.fluxDict:
            self.fluxDict[name] = 0.

    def setByModel(self, model, fluxVec):
        """ assign the fluxes given by fluxVec to the reactions in the given
            model - fluxVec[i] is assigned to the i-th reaction of the model;
            bounds are also set from the model

        Keyword arguments:

        model   -- a MetabolicModel
        fluxVec -- vector of fluxes, length must be number of reactions in the
                   model

        This function does not clear self.fluxDict beforehand.
        """
        if len(fluxVec) == 0:
            return  # Nothing to do (e.g. if FBA solution is empty)
        if len(model) != len(fluxVec):
            raise ValueError("Length of flux vector and number of reactions "
                             "in the model don't agree")
        for (name, flux) in zip(model.getReactionNames(), fluxVec):
            self.fluxDict[name] = flux
        lb, ub = model.getBounds()
        for name, bounds in zip(model.getReactionNames(), list(zip(lb, ub))):
            self.boundsDict[name] = bounds

    def setByDict(self, d, fluxVec):
        """ set the fluxes from the given vector with the order from the given
            dictionary

        Keyword arguments:

        d       -- dictionary {name : index} for ordering fluxVec
        fluxVec -- vector of fluxes, length must be number of reactions in the
                   model

        This function does not clear self.fluxDict beforehand.
        """
        if len(d) != len(fluxVec):
            raise ValueError(
                "Length of flux vector and dictionary don't agree")
        for name in d:
            self.fluxDict[name] = fluxVec[d[name]]

    def setByList(self, names, fluxVec):
        """ set the fluxes from the given vector with the order from the given
            dictionary

        Keyword arguments:

        d       -- dictionary {name : index} for ordering fluxVec
        fluxVec -- vector of fluxes, length must be number of reactions in the
                   model

        This function does not clear self.fluxDict beforehand.
        """
        if len(names) != len(fluxVec):
            raise ValueError("Length of name list and flux vector don't agree")
        self.fluxDict = dict(list(zip(names, fluxVec)))

    def getVecOrderedByModel(self, model):
        """ return a flux vector ordered like the reactions in the model

        Fluxes through reactions not present in the model are ignored.
        """
        fluxVec = []
        for name in model.getReactionNames():
            if name in self.fluxDict:
                fluxVec.append(self.fluxDict[name])
            else:
                fluxVec.append(None)  # Reaction has no flux value
        return fluxVec

    def getVecOrderedByDict(self, d):
        """ return a flux vector with the order from the given dictionary

        Keyword arguments:

        d       -- dictionary {name : index} for ordering fluxVec
        """
        fluxVec = [0.]*len(d)
        for name in d:
            if name in self.fluxDict:
                fluxVec[d[name]] = self.fluxDict[name]
            else:
                fluxVec[d[name]] = None  # Reaction has no flux value
        return fluxVec

    def hasSameReactions(self, reactions):
        """ verify that the MetabolicFlux object has the same reactions as the
            given 'reactions' object, which can be a MetabolicModel, a dict,
            a sequence type or a MetabolicFlux object
        """
        # Get set of reactions (to eliminate possible duplicates)
        if isinstance(reactions, MetabolicModel):
            reactionSet = set(reactions.getReactionNames())
        else:
            reactionSet = set(reactions)

        if len(self.fluxDict) != len(reactionSet):
            return False
        for rea in reactionSet:
            if rea not in self.fluxDict:
                return False
        return True

    def absDiff(self, other):
        """ return the absolute difference to the other MetabolicFlux as
            MetabolicFlux object; if a reaction is only contained in self or
            other, the flux is treated as if it were zero in the other object
        """
        diffDict = {}
        for rea in self.fluxDict:
            if rea in other.fluxDict:
                diffDict[rea] = abs(self.fluxDict[rea]-other.fluxDict[rea])
            else:
                diffDict[rea] = abs(self.fluxDict[rea])

        for rea in other.fluxDict:
            if rea not in self.fluxDict:
                diffDict[rea] = abs(other.fluxDict[rea])
        return MetabolicFlux(diffDict)

    def strictAbsDiff(self, other):
        """ return the absolute difference to the other MetabolicFlux as
            MetabolicFlux object - both must have exactly the same reactions
        """
        diffDict = {}
        if len(self.fluxDict) != len(other.fluxDict):
            raise ValueError("Lengths don't agree")
        for name in self.fluxDict:
            if name not in other.fluxDict:
                raise ValueError("Reaction '%s' not found in other "
                                 "MetabolicFlux object" % name)
            diffDict[name] = abs(self.fluxDict[name]-other.fluxDict[name])
        return MetabolicFlux(diffDict)

    def sqDiff(self, other):
        """ return the squared difference to the other MetabolicFlux as
            MetabolicFlux object; if a reaction is only contained in self or
            other, the flux is treated as if it were zero in the other object
        """
        diffDict = {}
        for rea in self.fluxDict:
            if rea in other.fluxDict:
                diff = self.fluxDict[rea]-other.fluxDict[rea]
            else:
                diff = self.fluxDict[rea]
            diffDict[rea] = diff*diff

        for rea in other.fluxDict:
            if rea not in self.fluxDict:
                diff = other.fluxDict[rea]
                diffDict[rea] = diff*diff
        return MetabolicFlux(diffDict)

    def strictSqDiff(self, other):
        """ return the squared difference to the other MetabolicFlux as
            MetabolicFlux object - both must have exactly the same reactions
        """
        diffDict = {}
        if len(self.fluxDict) != len(other.fluxDict):
            raise ValueError("Lengths don't agree")
        for name in self.fluxDict:
            if name not in other.fluxDict:
                raise ValueError("Reaction '%s' not found in other "
                                 "MetabolicFlux object" % name)
            diff = self.fluxDict[name]-other.fluxDict[name]
            diffDict[name] = diff*diff
        return MetabolicFlux(diffDict)

    def computeSplitRatios(self, metabolite, model, cutoff=0.,
                           cutoffIsAbsolute=False, listAll=False):
        """ compute split ratios for the fluxes entering and leaving the given
            metabolite node

        Keyword arguments:

        metabolite       -- name of the metabolite
        model            -- the underlying MetabolicModel
        cutoff           -- cutoff value below which fluxes are considered
                            absent, must be >= 0
        cutoffIsAbsolute -- True : cutoff is in absolute flux units,
                            False: cutoff is relative (between 0 and 1; default)
        listAll          -- if True, ignore cutoff and list all fluxes (even 0)

        Returns:
        outRatios, inRatios

        outRatios        -- dict { reaction : (ratio, flux) } for all outgoing
                            fluxes, i.e. the fluxes consuming the metabolite
        inRatios         -- dict { reaction : (ratio, flux) } for all incoming
                            fluxes, i.e. the fluxes producing the metabolite
        """
        posDict, negDict = {}, {}
        for rea in model:
            for coef, met in rea:
                if met == metabolite:
                    # Compute flux through metabolite from current reaction
                    metFlux = self.fluxDict[rea.name]*coef
                    if listAll:
                        # Zero fluxes are arbitrarily counted as outgoing
                        if metFlux > 0.:
                            addToDictPlus(posDict, rea.name, metFlux)
                        else:
                            addToDictPlus(negDict, rea.name, -metFlux)
                    if cutoffIsAbsolute:
                        if metFlux > cutoff:
                            addToDictPlus(posDict, rea.name, metFlux)
                        elif metFlux < -cutoff:
                            addToDictPlus(negDict, rea.name, -metFlux)
                    else:
                        # Fluxes of exactly zero can be discarded already here
                        if metFlux > 0.:
                            addToDictPlus(posDict, rea.name, metFlux)
                        elif metFlux < 0.:
                            addToDictPlus(negDict, rea.name, -metFlux)

        dontCheck = listAll or cutoffIsAbsolute

        posSum = sum(posDict.values())
        reactions = list(posDict.keys())
        for rea in reactions:
            value = posDict[rea]
            ratio = 0. if posSum == 0. else old_div(value, posSum)
            if dontCheck or ratio > cutoff:
                posDict[rea] = ratio, value
            else:
                del posDict[rea]

        negSum = sum(negDict.values())
        reactions = list(negDict.keys())
        for rea in reactions:
            value = negDict[rea]
            ratio = 0. if negSum == 0. else old_div(value, negSum)
            if dontCheck or ratio > cutoff:
                negDict[rea] = ratio, value
            else:
                del negDict[rea]

        return negDict, posDict

    def computeAllSplitRatios(self, model, cutoff=0., cutoffIsAbsolute=False,
                              listAll=False):
        """ compute metabolite fluxes and split ratios for the fluxes entering
            and leaving each metabolite node in the given model

        Keyword arguments:

        model            -- the underlying MetabolicModel
        cutoff           -- cutoff value below which fluxes are considered
                            absent, must be >= 0
        cutoffIsAbsolute -- True : cutoff is in absolute flux units,
                            False: cutoff is relative (between 0 and 1; default)
        listAll          -- if True, ignore cutoff and list all fluxes (even 0)

        Returns:
        splitRatios      -- dict { metabolite : (outRatios, inRatios) } with
                            outRatios and inRatios as returned by
                            computeSplitRatios()
        """
        splitRatios = dict((met, ({}, {}))
                           for met in model.getMetaboliteNames())
        for rea in model:
            for coef, met in rea:
                # Compute flux through metabolite 'met' from reaction 'rea'
                metFlux = self.fluxDict[rea.name]*coef
                if listAll:
                    # Zero fluxes are arbitrarily counted as outgoing
                    if metFlux > 0.:
                        addToDictPlus(splitRatios[met][1], rea.name, metFlux)
                    else:
                        addToDictPlus(splitRatios[met][0], rea.name, -metFlux)
                elif cutoffIsAbsolute:
                    if metFlux > cutoff:
                        addToDictPlus(splitRatios[met][1], rea.name, metFlux)
                    elif metFlux < -cutoff:
                        addToDictPlus(splitRatios[met][0], rea.name, -metFlux)
                else:
                    # Fluxes of exactly zero can be discarded already here
                    if metFlux > 0.:
                        addToDictPlus(splitRatios[met][1], rea.name, metFlux)
                    elif metFlux < 0.:
                        addToDictPlus(splitRatios[met][0], rea.name, -metFlux)

        dontCheck = listAll or cutoffIsAbsolute

        for met in splitRatios:
            negDict, posDict = splitRatios[met]
            posSum = sum(posDict.values())
            reactions = list(posDict.keys())
            for rea in reactions:
                value = posDict[rea]
                ratio = 0. if posSum == 0. else old_div(value, posSum)
                if dontCheck or ratio > cutoff:
                    posDict[rea] = ratio, value
                else:
                    del posDict[rea]

            negSum = sum(negDict.values())
            reactions = list(negDict.keys())
            for rea in reactions:
                value = negDict[rea]
                ratio = 0. if negSum == 0. else old_div(value, negSum)
                if dontCheck or ratio > cutoff:
                    negDict[rea] = ratio, value
                else:
                    del negDict[rea]

        return splitRatios

    # ------ File I/O ----------------------------------------------------------

    def writeToFile(self, filename):
        """ write to the given file
        """
        with open(filename, 'w') as f:
            if filename.split(".")[-1] == "csv":
                self.writeToCSVHandle(f)
            else:
                self.writeToFileHandle(f)

    def writeToCSVHandle(self, f):
        """ write to the CSV file given by file handle f (must be open for writing)
        """
        if not self.fluxDict:
            return  # Nothing to do

        csvwriter = csv.writer(
            f, delimiter=";", quotechar='"', quoting=csv.QUOTE_MINIMAL)
        # Write Header to CSV file
        csvwriter.writerow(["NAME", "FLUX", "LB", "UB"])
        result = []
        # Construct output table and get widths of table columns (length of
        # longest occurring string)
        for rea in self.fluxDict:
            try:
                lb, ub = self.boundsDict[rea]
            except KeyError:
                lb, ub = -inf, inf
            # update python3
            result.append([rea.name if type(rea) != str else rea,
                           self.fluxDict[rea], lb, ub])

        for row in sorted(result):
            csvwriter.writerow(row)

    def writeToFileHandle(self, f):
        """ write to the file given by file handle f (must be open for writing)
        """
        if not self.fluxDict:
            return  # Nothing to do

        # Construct output table and get widths of table columns (length of
        # longest occurring string)
        maxlenName, maxlenFlux, maxlenLb, maxlenUb = 0, 0, 0, 0
        outputTable = {}  # dict {name : (flux, lb, ub)} - all as strings
        for rea in self.fluxDict:
            try:
                lb, ub = self.boundsDict[rea]
            except KeyError:
                lb, ub = -inf, inf

            # Pad non-negative numbers with space
            outputTable[rea] = tuple(map(padNumber, list(map(repr,
                                                             (self.fluxDict[rea], lb, ub)))))
            lenName, lenFlux, lenLb, lenUb = list(
                map(len, (rea,)+outputTable[rea]))

            if lenName > maxlenName:
                maxlenName = lenName
            if lenFlux > maxlenFlux:
                maxlenFlux = lenFlux
            if lenLb > maxlenLb:
                maxlenLb = lenLb
            if lenUb > maxlenUb:
                maxlenUb = lenUb

        # Write table head
        f.write("NAME".ljust(maxlenName)+"   "+"FLUX".center(maxlenFlux)+" " +
                "LB  ".rjust(maxlenLb)+" "+"UB  ".rjust(maxlenUb)+"\n")
        # Write result for the flux through every reaction
        for rea in sorted(self.fluxDict):
            flux, lb, ub = outputTable[rea]
            # update python3
            if type(rea) != str:
                rea = rea.name
            f.write(rea.ljust(maxlenName) + " : " + flux.rjust(maxlenFlux) +
                    " " + lb.rjust(maxlenLb) + " " + ub.rjust(maxlenUb) + "\n")

    def readFromFile(self, filename):
        """ parse the given file
        """
        with open(filename, 'r') as f:
            if filename.split(".")[-1] == "csv":
                self.readFromCSVHandle(f)
            else:
                self.readFromFileHandle(f)

    def readFromCSVHandle(self, f):
        """ parse the file given by file handle f (must be open for reading)
        """
        csvreader = csv.reader(f, delimiter=";", quotechar='"')

        reactions = set()
        line_no = 0
        for line in csvreader:
            line_no += 1
            if line[0] == "NAME":
                continue

            try:
                self.fluxDict[line[0]] = float(line[1])
                self.boundsDict[line[0]] = line[2], line[3]
            except IndexError:
                raise SyntaxError("Syntax error in line %u:\nLine must "
                                  "contain exactly four values (name, flux, "
                                  "lb, ub)." % line_no)
            except ValueError:
                raise SyntaxError("Syntax error in line %u: Invalid "
                                  "floating point value." % line_no)

    def readFromFileHandle(self, f):
        """ parse the file given by file handle f (must be open for reading)
        """
        reactions = set()
        line_no = 0
        for line in f:
            line_no += 1
            if (line.lstrip().upper().startswith("NAME") or line == "" or
                    line.isspace()):
                continue

            try:
                rea, values_str = list(map(str.rstrip, line.split(":")))
            except ValueError:
                raise SyntaxError("Syntax error in line %u:\nLine must "
                                  "contain exactly one colon (':')." %
                                  line_no)
            if rea in reactions:
                raise SyntaxError("Syntax error in line %u: Duplicate "
                                  "reaction." % line_no)
            reactions.add(rea)

            try:
                values = list(map(float, values_str.split(None, 3)[:3]))
            except ValueError:
                raise SyntaxError("Syntax error in line %u: Invalid "
                                  "floating point value." % line_no)
            try:
                self.fluxDict[rea] = values[0]
                self.boundsDict[rea] = values[1], values[2]
            except IndexError:
                raise SyntaxError("Syntax error in line %u:\nLine must "
                                  "contain exactly three values (flux, lb, "
                                  "ub)." % line_no)
