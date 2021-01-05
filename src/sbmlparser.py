#!/usr/bin/env python
"""
This module defines a parser for SBML files. It is also a script for converting
SBML files to metano files. This parser ignores units. The most recent version
supported is SBML Level 2 Version 4, with the Flux Balance Constraints
extension: http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_(flux).


This file is part of metano.
Copyright (C) 2010-2019 Alexander Riemer, Julia Helmecke
Braunschweig University of Technology,
Dept. of Bioinformatics and Biochemistry

Alexander Krause contributed to this file.

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

from builtins import zip
from builtins import object
from metano.defines import (Verbosity, ReaFileStruc, makeCoefNameTerm, isZeroVector,
                     COPYRIGHT_VERSION_STRING)
from metano.sbmlstructure import SbmlDef, getFirstElementChild
from metano.fba import OptionParser
from metano.metabolicmodel import Reactant, Reaction, MetabolicModel
from metano.metabolicflux import MetabolicFlux
import xml.dom.minidom as dom
from xml.parsers.expat import ExpatError
from numpy import inf
import os, re


class SbmlParser(object):
    """ Class for parsing SBML files

    Member variables:

    modelName          -- name of SBML model
    listOfReactions    -- list of Reaction objects (class Reaction defined in
                          metabolicmodel.py)
    listOfFluxValues   -- list of fluxes corresponding to the reactions in
                          listOfReactions (left empty, if no fluxes are found)
    listOfObjCoef      -- list of OBJECTIVE_COEFFICIENT values for each reaction
    setOfReactions     -- IDs of reactions already read
    dictOfCompartments -- dict { ID : (name, outside) } for each compartment
    dictOfSpecies      -- stores tuple (name, compartment, boundaryCondition)
                          for each species ID
    dictOfReactions    -- dict { ID : name } for each reaction
    infLevel           -- threshold above which a bound is considered infinite
                          (default: 1000)
    messages           -- error, warning, info, and debug messages generated
                          during parsing
    """

    _decodeBool = {"false": 0, "False": 0, "0": 0, "true": 1, "True": 1, "1": 1}


    def __init__(self, infLevel=1000.):
        """ initialize class SbmlParser
        """
        self.infLevel = infLevel

    def clear(self):
        """ clear all internal variables, warnings and info messages
        """
        self.messages = []
        self.modelName = ""
        self.listOfReactions = []
        self.setOfReactions = set()
        self.listOfFluxValues = []
        self.listOfObjCoef = []
        self.dictOfCompartments = {}
        self.dictOfSpecies = {}
        self.dictOfReactions = {}
        self._tmpReaction = CacheReaction()

    def getMetabolicModelComponents(self):
        """ return structures needed to build a MetabolicModel
        """
        # Return only non-boundary species, with 'isBoundary' field removed
        metaboliteLabels = dict((ID, (name, compartment)) for
                                (ID, (name, compartment, isBoundary)) in
                                list(self.dictOfSpecies.items()) if not isBoundary)

        return self.listOfReactions, metaboliteLabels, self.dictOfReactions

    def getFluxValues(self, reaNameDict=None):
        """ return the flux distribution read from the SBML file as dict
            { reaction ID/name : flux value }

        Keyword arguments:

        reaNameDict -- dict {reaction ID : name } (optional; if present,
                                                   names must be unique)
        """
        if reaNameDict:
            return dict(list(zip([reaNameDict[x.name] for x in self.listOfReactions],
                            self.listOfFluxValues)))
        else:
            return dict(list(zip([x.name for x in self.listOfReactions],
                            self.listOfFluxValues)))

    def getObjCoefficients(self, reaNameDict=None):
        """ return the coefficient vector of the linear objective function read
            from the SBML file as dict { reaction ID/name : coefficient }

        Keyword arguments:

        reaNameDict -- dict {reaction ID : name } (optional; if present,
                                                   names must be unique)
        """
        if reaNameDict:
            return dict(list(zip([reaNameDict[x.name] for x in self.listOfReactions],
                            self.listOfObjCoef)))
        else:
            return dict(list(zip([x.name for x in self.listOfReactions],
                            self.listOfObjCoef)))

    def getMessages(self, level=Verbosity.INFO):
        """ return a list of all messages at or above the given level of
            severity

        Levels are defined in defines.py.
        """
        return [x[1] for x in self.messages if x[0] <= level]


    # Shortcuts for generating common warnings, error and debug messages

    def _errorEmptyCollection(self, collectionName):
        self.messages.append((Verbosity.ERROR,
                              "Error: %s is empty." % collectionName))
    def _errorNoId(self, nodeType):
        self.messages.append((Verbosity.ERROR,
            "Error: Found %s with no ID or name." % nodeType))
    def _errorDuplicateId(self, ID, collectionName):
        self.messages.append((Verbosity.ERROR,
            "Error: Found duplicate ID (%s) in %s." % (ID, collectionName)))
    def _errorIllegalFloat(self, node, upperType, upperId):
        self.messages.append((Verbosity.ERROR, "Error: Illegal floating point "
                              "value for %s in %s '%s'." % (node, upperType,
                                                            upperId)))
    def _errorMissingValue(self, node, upperType, upperId):
        self.messages.append((Verbosity.ERROR, "Error in %s '%s': Parameter "
            "%s has no %s." % (upperType, upperId, node, SbmlDef.VALUE)))

    def _warnDuplicateNode(self, node, upperType="", upperId=""):
        msg = "Warning: Duplicate %s found" % node
        if upperId:
            msg += " in %s '%s'" % (upperType, upperId)
        msg += ". Only the first is used."
        self.messages.append((Verbosity.WARNING, msg))
    def _warnIllegalValueUsingDefault(self, nodeType, ID, attribute, default):
        self.messages.append((Verbosity.WARNING,
            "Warning: Unknown value encountered for attribute \"%s\" in %s "
            "'%s'. Set to default: %s." % (attribute, nodeType, ID, default)))

    def _debugAttributeValue(self, attribute, value):
        self.messages.append((Verbosity.DEBUG,
                              "Debug: %s = %s." % (attribute, value)))
    def _debugUsingIdAsName(self, nodeType, ID, nameStr):
        self.messages.append((Verbosity.DEBUG, "Debug: %s '%s' is missing "
                              "attribute \"%s\". Using ID instead." %
                              (nodeType.capitalize(), ID, nameStr)))
    def _debugUsingNameAsId(self, nodeType, name):
        self.messages.append((Verbosity.DEBUG, "Debug: %s named '%s' has no ID."
                              " Using name instead." %
                              (nodeType.capitalize(), name)))
    def _debugIgnoringNode(self, node, upperId="", upperType=""):
        msg = "Debug: Ignoring SBML node \"%s\"" % node
        if upperId:
            if upperType:
                msg += " under %s '%s'" % (upperType, upperId)
            else:
                msg += " in %s" % upperId
        msg += "."
        self.messages.append((Verbosity.DEBUG, msg))


    # Main methods

    def readFile(self, filename, boundaryRegex=""):
        """ read the given SBML file into member structures

        Keyword arguments:

        filename      -- name of the SBML file to be parsed
        boundaryRegex -- if given, species with ID matching this regular
                         expression are treated as additional boundary species

        Returns True if successful, else False
        """
        self.clear()
        isModelFound = False
        fc = dom.parse(filename)

        # Read model node first
        # There are ELEMENT_NODES with subtrees and TEXT_NODES - only the
        # ELEMENT_NODES are processed
        for node in getFirstElementChild(fc).childNodes:
            if node.nodeName.lower() == SbmlDef.MODEL.lower():
                isModelFound = True

                # Read model name (optional)
                if node.getAttribute(SbmlDef.NAME) != "":
                    self.modelName = node.getAttribute(SbmlDef.NAME)
                    msg = ("Info: The metabolic model has the name '%s'."
                           % self.modelName)
                    self.messages.append((Verbosity.INFO, msg))
                else:
                    msg = "Info: The metabolic model has no name."
                    self.messages.append((Verbosity.INFO, msg))

                # Read the model elements
                success = self._readModelNodes(node)
                break  # an SBML file may contain only one model

        if not isModelFound:
            success = False
            msg = "Error: Unable to read model."
            self.messages.append((Verbosity.ERROR, msg))

        if success:
            nBoundary = len([isBoundary for (_, _, isBoundary) in
                             list(self.dictOfSpecies.values()) if isBoundary])
            msg = ("Info: The SBML model has %d reactions and %d species" %
                   (len(self.listOfReactions), len(self.dictOfSpecies)))
            if nBoundary:
                msg += " (including %d boundary species)" % nBoundary
            msg +="."
            self.messages.append((Verbosity.INFO, msg))

            if boundaryRegex:
                try:
                    pattern = re.compile(boundaryRegex)
                except re.error as strerror:
                    msg = ("Warning: Error in regular expression for boundary "
                           "metabolites: %s. Unable to evaluate." %
                           strerror)
                    self.messages.append((Verbosity.WARNING, msg))
                    boundaryRegex = ""
                    pattern = None

                if pattern:
                    # Delete all metabolites matching the given reg-ex
                    for reaction in self.listOfReactions:
                        for reactant in reaction:
                            metName = reactant.name
                            if pattern.search(metName):
                                reaction.removeReactant(metName)
                                msg = ("Debug: Deleted boundary metabolite '%s'"
                                       " from reaction '%s'." % (metName,
                                                                 reaction.name))
                                self.messages.append((Verbosity.DEBUG, msg))

        return success


    def fixNames(self):
        """ modify reaction, species, and compartment IDs and names by replacing
            all characters not allowed in metano reaction files with underscores
            and converting all strings to type str
        """
        illegC = re.compile("\\s|%s|%s" % (ReaFileStruc.delimiter,
                                           ReaFileStruc.commentSign))
        # Reaction IDs and species IDs in reactions
        for reaction in self.listOfReactions:
            reaction.name = re.sub(illegC, "_", reaction.name).encode('ascii',
                                                                      'replace')
            for reactant in reaction:
                reactant.name = re.sub(illegC, "_", reactant.name).encode(
                                                   'ascii', 'replace')

        # Reaction names
        strDictOfReactions = {}
        for rea in self.dictOfReactions:
            strDictOfReactions[re.sub(illegC, "_", rea).encode('ascii',
                'replace')] = re.sub(illegC, "_",
                self.dictOfReactions[rea]).encode('ascii', 'replace')

        self.dictOfReactions = strDictOfReactions

        strDictOfSpecies = {}
        for met in self.dictOfSpecies:
            # Species names and compartment IDs in species
            name, compartment, isBoundary = self.dictOfSpecies[met]
            strDictOfSpecies[re.sub(illegC, "_", met).encode('ascii',
                'replace')] = (re.sub(illegC, "_", name).encode('ascii',
                'replace'), re.sub(illegC, "_", compartment).encode('ascii',
                'replace'), isBoundary)

        self.dictOfSpecies = strDictOfSpecies

        strDictOfCompartments = {}
        for cId in self.dictOfCompartments:
            # Name of compartment and ID of outside compartment
            name, outside = self.dictOfCompartments[cId]
            strDictOfCompartments[re.sub(illegC, "_", cId).encode('ascii',
                'replace')] = (re.sub(illegC, "_", name).encode('ascii',
                'replace'), re.sub(illegC, "_", outside).encode('ascii',
                                                                'replace'))

        self.dictOfCompartments = strDictOfCompartments


    def buildModel(self, usenames=False):
        """ return a MetabolicModel constructed from the data read from the SBML
            file

        Keyword arguments:

        usenames -- if True, write reaction and species names as labels to the
                    model
        """
        reactions, metaboliteLabels, reactionLabels = \
            self.getMetabolicModelComponents()

        model = MetabolicModel()
        model.addReactions(reactions)

        # Set labels built from name and compartment if names are to be used
        # instead of IDs
        if usenames:
            metLabelDict = dict((k, "%s[%s]" % x)
                                for (k, x) in list(metaboliteLabels.items()))
            model.setMetaboliteLabels(metLabelDict)
            model.setReactionLabels(reactionLabels)
            # Make sure that all reaction and metabolite labels are unique
            model.makeLabelsUnique()

        return model


    # Auxiliary methods called by readFile() to parse subtrees

    def _readModelNodes(self, modelNode):
        """ read the subtree under the model node

        Returns True if successful, else False
        """
        isListOfCompartmentsFound = False
        isListOfCompartmentsRead = False
        isListOfSpeciesFound = False
        isListOfSpeciesRead = False
        isListOfReactionsFound = False
        isListOfReactionsRead = False

        # Notes:
        # 1) If the model contains more than one listOfCompartments,
        #    listOfSpecies, or listOfReactions, only the first occurrence of
        #    each is processed.
        # 2) Order is important in SBML: listOfSpecies cannot occur before
        #    listOfCompartments, and listOfReactions cannot occur before
        #    listOfSpecies

        for node in modelNode.childNodes:
            nodeName = node.nodeName.lower()

            if nodeName == SbmlDef.LISTOFCOMPARTMENTS.lower():

                if isListOfCompartmentsFound:
                    self._warnDuplicateNode(SbmlDef.LISTOFCOMPARTMENTS)
                else:
                    isListOfCompartmentsFound = True
                    isListOfCompartmentsRead = \
                        self._readListOfCompartments(node)

            elif nodeName == SbmlDef.LISTOFSPECIES.lower():

                if isListOfSpeciesFound:
                    self._warnDuplicateNode(SbmlDef.LISTOFSPECIES)

                elif not isListOfCompartmentsRead:
                    msg = ("Error: Encountered %s, but need %s first. "
                           "Unable to process." % (SbmlDef.LISTOFSPECIES,
                                                   SbmlDef.LISTOFCOMPARTMENTS))
                    self.messages.append((Verbosity.ERROR, msg))
                    isListOfSpeciesFound = True
                else:
                    isListOfSpeciesFound = True
                    isListOfSpeciesRead = self._readListOfSpecies(node)

            elif nodeName == SbmlDef.LISTOFREACTIONS.lower():

                if isListOfReactionsFound:
                    self._warnDuplicateNode(SbmlDef.LISTOFREACTIONS)

                elif not isListOfSpeciesRead:
                    msg = ("Error: Encountered %s, but need %s first. "
                           "Unable to process." % (SbmlDef.LISTOFREACTIONS,
                                                   SbmlDef.LISTOFSPECIES))
                    self.messages.append((Verbosity.ERROR, msg))
                    isListOfReactionsFound = True
                else:
                    isListOfReactionsFound = True
                    isListOfReactionsRead = self._readListOfReactions(node)

            elif node.nodeType == dom.Node.ELEMENT_NODE and nodeName != "":
                # Anything that is not a listOfCompartments, listOfSpecies, or
                # listOfReactions is ignored
                self._debugIgnoringNode(node.nodeName)

        if (isListOfSpeciesRead and isListOfReactionsRead and
            isListOfCompartmentsRead):
            return True

        if not isListOfCompartmentsFound:
            self.messages.append((Verbosity.ERROR,
                "Error: %s not found." % SbmlDef.LISTOFCOMPARTMENTS))
        elif not isListOfSpeciesFound:
            self.messages.append((Verbosity.ERROR,
                "Error: %s not found." % SbmlDef.LISTOFSPECIES))
        elif not isListOfReactionsFound:
            self.messages.append((Verbosity.ERROR,
                "Error: %s not found." % SbmlDef.LISTOFREACTIONS))

        return False


    def _readListOfCompartments(self, compartmentsNode):
        """ read the subtree under the listOfCompartments node

        Returns True if successful, else False
        """
        success = True
        for node in compartmentsNode.childNodes:

            if node.nodeName.lower() == SbmlDef.COMPARTMENT.lower():
                compartmentOutside = ""
                compartmentName = ""

                if node.getAttribute(SbmlDef.NAME) != "":
                    compartmentName = node.getAttribute(SbmlDef.NAME)

                if node.getAttribute(SbmlDef.ID) != "":
                    compartmentId = node.getAttribute(SbmlDef.ID)
                elif compartmentName != "":
                    compartmentId = compartmentName
                    self._debugUsingNameAsId("compartment", compartmentName)
                else:
                    self._errorNoId("compartment")
                    success = False
                    continue

                if node.getAttribute(SbmlDef.COMPARTMENTOUTSIDE) != "":
                    compartmentOutside = \
                        node.getAttribute(SbmlDef.COMPARTMENTOUTSIDE)

                # ID must be unique
                if compartmentId in self.dictOfCompartments:
                    self._errorDuplicateId(compartmentId,
                                           SbmlDef.LISTOFCOMPARTMENTS)
                    success = False
                    continue

                self.dictOfCompartments[compartmentId] = (compartmentName,
                                                          compartmentOutside)

                msg = "Info: Found compartment '%s'" % compartmentId
                if compartmentName:
                    msg += " (%s)" % compartmentName
                if compartmentOutside:
                    msg += " with %s='%s'" % (SbmlDef.COMPARTMENTOUTSIDE,
                                              compartmentOutside)
                msg += "."
                self.messages.append((Verbosity.INFO, msg))

            elif node.nodeType == dom.Node.ELEMENT_NODE and node.nodeName != "":
                # Anything that is not a compartment is ignored
                self._debugIgnoringNode(node.nodeName,
                                        SbmlDef.LISTOFCOMPARTMENTS)

        if not self.dictOfCompartments:
            self._errorEmptyCollection(SbmlDef.LISTOFCOMPARTMENTS)
            success = False

        return success


    def _readListOfSpecies(self, listOfSpeciesNode):
        """ read the subtree under the listOfSpecies node

        Returns True if successful, else False
        """
        success = True
        for node in listOfSpeciesNode.childNodes:

            if node.nodeName.lower() == SbmlDef.SPECIES.lower():
                speciesName = ""
                speciesCompartment = ""
                speciesBoundaryCondition = 0

                if node.getAttribute(SbmlDef.NAME) != "":
                    speciesName = node.getAttribute(SbmlDef.NAME)

                # Species must always have ID/name and compartment
                if node.getAttribute(SbmlDef.ID) != "":
                    speciesId = node.getAttribute(SbmlDef.ID)
                elif speciesName != "":
                    speciesId = speciesName
                    self._debugUsingNameAsId("species", speciesName)
                else:
                    self._errorNoId("species")
                    success = False
                    continue

                # Species ID must be unique
                if speciesId in self.dictOfSpecies:
                    self._errorDuplicateId(speciesId, SbmlDef.LISTOFSPECIES)
                    success = False
                    continue

                if node.getAttribute(SbmlDef.COMPARTMENT) != "":
                    speciesCompartment = node.getAttribute(SbmlDef.COMPARTMENT)
                else:
                    msg = ("Error: Species '%s' is missing attribute \"%s\"."
                           % (speciesId, SbmlDef.COMPARTMENT))
                    self.messages.append((Verbosity.ERROR, msg))
                    success = False
                    continue

                msg = "Debug: Found species '%s'" % speciesId
                if speciesName:
                    msg += " (%s)" % speciesName
                msg += " in %s '%s'." % (SbmlDef.COMPARTMENT,
                                         speciesCompartment)
                self.messages.append((Verbosity.DEBUG, msg))

                if not speciesName:
                    speciesName = speciesId
                    self._debugUsingIdAsName("species", speciesId,
                                             SbmlDef.NAME)

                # Process Boolean attribute 'boundaryCondition'
                bcStr = node.getAttribute(SbmlDef.BOUNDARYCONDITION)
                if bcStr != "":
                    if bcStr in self._decodeBool:
                        speciesBoundaryCondition = self._decodeBool[bcStr]
                        self._debugAttributeValue(SbmlDef.BOUNDARYCONDITION,
                                                  speciesBoundaryCondition)
                    else:
                        self._warnIllegalValueUsingDefault("species", speciesId,
                            SbmlDef.BOUNDARYCONDITION, False)

                # Is compartment known?
                if speciesCompartment not in self.dictOfCompartments:
                    msg = ("Error: Undefined %s (%s) found for species '%s.'" %
                           (SbmlDef.COMPARTMENT, speciesCompartment, speciesId))
                    self.messages.append((Verbosity.ERROR, msg))
                    success = False
                    continue

                self.dictOfSpecies[speciesId] = (speciesName,
                               speciesCompartment, speciesBoundaryCondition)

            elif node.nodeType == dom.Node.ELEMENT_NODE and node.nodeName != "":
                # Anything that is not a species is ignored
                self._debugIgnoringNode(node.nodeName, SbmlDef.LISTOFSPECIES)

        if not self.dictOfSpecies:
            self._errorEmptyCollection(SbmlDef.LISTOFSPECIES)
            success = False

        return success


    def _readListOfReactions(self, listOfReactionsNode):
        """ read subtree under the listOfReactions node

        Returns True if successful, else False
        """
        success = True
        listIsEmpty = True
        for node in listOfReactionsNode.childNodes:
            if listIsEmpty:
                listIsEmpty = False
            if node.nodeName.lower() == SbmlDef.REACTION.lower():
                if not self._readReaction(node):
                    success = False

            elif node.nodeType == dom.Node.ELEMENT_NODE and node.nodeName != "":
                # Anything that is not a reaction is ignored
                self._debugIgnoringNode(node.nodeName, SbmlDef.LISTOFREACTIONS)

        if listIsEmpty:
            self._errorEmptyCollection(SbmlDef.LISTOFREACTIONS)
            success = False

        return success


    def _readReaction(self, reactionNode):
        """ process reaction node

        Returns True if successful, else False
        """
        success = True
        reactionId = ""
        reactionName = ""
        isReversible = 1
        isListOfReactantsFound = False
        isListOfProductsFound = False
        isKineticLawFound = False

        # Read reaction ID or generate if not present
        if reactionNode.getAttribute(SbmlDef.ID) != "":
            reactionId = reactionNode.getAttribute(SbmlDef.ID)
            msg = "Debug: Reaction ID is '%s'." % reactionId
            self.messages.append((Verbosity.DEBUG, msg))
        else:
            reactionId = "rea%03d" % (len(self.listOfReactions)+1)
            msg = ("Warning: Found reaction without attribute \"%s\". Generated"
                   " ID: '%s'." % (SbmlDef.ID, reactionId))
            self.messages.append((Verbosity.WARNING, msg))

        # ID must be unique
        if reactionId in self.setOfReactions:
            self._errorDuplicateId(reactionId, SbmlDef.LISTOFREACTIONS)
            success = False
        else:
            self.setOfReactions.add(reactionId)

        # Read reaction name (optional)
        if reactionNode.getAttribute(SbmlDef.NAME) != "":
            reactionName = reactionNode.getAttribute(SbmlDef.NAME)
            msg = "Debug: Reaction name is '%s'." % reactionName
            self.messages.append((Verbosity.DEBUG, msg))
        else:
            reactionName = reactionId
            self._debugUsingIdAsName("reaction", reactionId, SbmlDef.NAME)

        # Process Boolean attribute 'reversible'
        reversibleStr = reactionNode.getAttribute(SbmlDef.REVERSIBLE)
        if reversibleStr != "":
            if reversibleStr in self._decodeBool:
                isReversible = self._decodeBool[reversibleStr]
                self._debugAttributeValue(SbmlDef.REVERSIBLE, isReversible)
            else:
                self._warnIllegalValueUsingDefault("reaction", reactionId,
                                                   SbmlDef.REVERSIBLE, True)
        else:
            msg = ("Debug: Attribute \"%s\" in reaction '%s' is empty or "
                   "missing. Set to default: True." % (SbmlDef.REVERSIBLE,
                                                       reactionId))
            self.messages.append((Verbosity.DEBUG, msg))

        self._tmpReaction.reactionId = reactionId
        self._tmpReaction.reactionName = reactionName
        self._tmpReaction.reversible = isReversible

        # Process child nodes (listOfReactants, listOfProducts, kineticLaw,
        # notes) - data is exchanged via self._tmpReaction
        for node in reactionNode.childNodes:
            nodeName = node.nodeName.lower()

            if nodeName == SbmlDef.LISTOFREACTANTS.lower():

                if isListOfReactantsFound:
                    self._warnDuplicateNode(SbmlDef.LISTOFREACTANTS, "reaction",
                                            reactionId)
                else:
                    isListOfReactantsFound = True
                    if not self._readListOfReactantsProducts(node, "reactant"):
                        success = False

            elif nodeName == SbmlDef.LISTOFPRODUCTS.lower():

                if isListOfProductsFound:
                    self._warnDuplicateNode(SbmlDef.LISTOFPRODUCTS, "reaction",
                                            reactionId)
                else:
                    isListOfProductsFound = True
                    if not self._readListOfReactantsProducts(node, "product"):
                        success = False

            elif nodeName == SbmlDef.KINETICLAW.lower():

                if isKineticLawFound:
                    self._warnDuplicateNode(SbmlDef.KINETICLAW, "reaction",
                                            reactionId)
                else:
                    isKineticLawFound = True
                    if not self._readKineticLaw(node):
                        success = False

            elif nodeName == SbmlDef.NOTES.lower():
                self._tmpReaction.notes = self._extractTextFromSubtree(node)

            elif node.nodeType == dom.Node.ELEMENT_NODE and nodeName != "":
                self._debugIgnoringNode(node.nodeName, reactionId, "reaction")

        isListOfReactantsFound = (isListOfReactantsFound and
                                  self._tmpReaction.reactants)
        isListOfProductsFound = (isListOfProductsFound and
                                 self._tmpReaction.products)

        # Check for completeness
        if not isListOfReactantsFound:
            if not isListOfProductsFound:
                msg = ("Error: Reaction '%s' has neither reactants nor "
                       "products." % reactionId)
                self.messages.append((Verbosity.ERROR, msg))
                success = False
            else:
                msg = ("Debug: Reaction '%s' has no reactants (left-hand side is"
                       " empty)." % reactionId)
                self.messages.append((Verbosity.DEBUG, msg))
        elif not isListOfProductsFound:
            msg = ("Debug: Reaction '%s' has no products (right-hand side is "
                   "empty)." % reactionId)
            self.messages.append((Verbosity.DEBUG, msg))

        if not isKineticLawFound:
            msg = ("Warning: Reaction '%s' has no parameters (bounds)." %
                   reactionId)
            self.messages.append((Verbosity.WARNING, msg))

        # If processing of child nodes was successful, create Reaction object
        if success:
            comment = self._tmpReaction.notes.strip()
            rea = Reaction(reactionId, self._tmpReaction.reactants,
                           products=self._tmpReaction.products, comment=comment,
                           lb=self._tmpReaction.lb, ub=self._tmpReaction.ub,
                           direction=1-self._tmpReaction.reversible)

            self.listOfReactions.append(rea)
            self.listOfFluxValues.append(self._tmpReaction.fluxVal)
            self.listOfObjCoef.append(self._tmpReaction.objCoef)
            self.dictOfReactions[reactionId] = reactionName

        self._tmpReaction.clear()
        return success


    def _readListOfReactantsProducts(self, listNode, typ):
        """ process listOfReactants or listOfProducts node

        typ is either 'reactant' or 'product'.

        Returns True if successful, else False
        """
        success = True
        reactionId = self._tmpReaction.reactionId
        reactantSet = set()

        for node in listNode.childNodes:
            if node.nodeName.lower() == SbmlDef.SPECIESREFERENCE.lower():
                speciesId = ""
                coef = 0.

                # Read relevant attributes (ID and stoichiometric coefficient)
                if not node.getAttribute(SbmlDef.SPECIES) == "":
                    speciesId = node.getAttribute(SbmlDef.SPECIES)
                else:
                    msg = ("Error: %s without \"%s\" attribute in reaction "
                           "'%s'." % (typ.capitalize(), SbmlDef.SPECIES,
                                      reactionId))
                    self.messages.append((Verbosity.ERROR, msg))
                    success = False
                    continue

                # Species ID must refer to an entry in listOfSpecies
                if speciesId not in self.dictOfSpecies:
                    msg = ("Error: %s '%s' encountered in reaction '%s' is not "
                           "present in %s." % (typ.capitalize(), speciesId,
                                             reactionId, SbmlDef.LISTOFSPECIES))
                    self.messages.append((Verbosity.ERROR, msg))
                    success = False
                    continue

                if node.getAttribute(SbmlDef.STOICHIOMETRY) != "":
                    try:
                        coef = float(node.getAttribute(SbmlDef.STOICHIOMETRY))
                    except ValueError:
                        msg = ("Error: %s '%s' in reaction '%s' has illegal "
                               "\"%s\" value (expected a float)." %
                               (typ.capitalize(), speciesId, reactionId,
                                SbmlDef.STOICHIOMETRY))
                        self.messages.append((Verbosity.ERROR, msg))
                        success = False
                    if coef < 0.:
                        msg = ("Error: Negative stoichiometric coefficient "
                               "encountered for %s '%s' in reaction '%s'." %
                               (typ, speciesId, reactionId))
                        self.messages.append((Verbosity.ERROR, msg))
                        success = False
                else:
                    coef = 1.
                    msg = ("Debug: Attribute \"%s\" of %s '%s' in reaction "
                           "'%s' is empty or missing. Set to default: 1." %
                           (SbmlDef.STOICHIOMETRY, typ, speciesId, reactionId))
                    self.messages.append((Verbosity.DEBUG, msg))

                # Skip reactant/product if it is a boundary species
                name, _, isBoundary = self.dictOfSpecies[speciesId]
                if isBoundary:
                    msg = ("Debug: Skipping boundary species '%s' (%s)." %
                           (speciesId, name))
                    self.messages.append((Verbosity.DEBUG, msg))
                else:
                    reactantSet.add(speciesId)
                    if success:
                        if typ == "reactant":
                            self._tmpReaction.addReactant(Reactant(speciesId,
                                                                   coef))
                        elif typ == "product":
                            self._tmpReaction.addProduct(Reactant(speciesId,
                                                                  coef))
                        msg = ("Debug: Added %s '%s' with coefficient %r." %
                               (typ, speciesId, coef))
                        self.messages.append((Verbosity.DEBUG, msg))

            elif node.nodeType == dom.Node.ELEMENT_NODE and node.nodeName != "":
                # Anything that is not a speciesReference is ignored
                if typ == "product":
                    collectionName = SbmlDef.LISTOFPRODUCTS
                else:
                    collectionName = SbmlDef.LISTOFREACTANTS
                self._debugIgnoringNode(node.nodeName, "%s of reaction '%s'" %
                                        collectionName, reactionId)

        return success


    @staticmethod
    def _extractTextFromSubtree(root):
        text = ""
        for node in root.childNodes:
            if node.nodeType == dom.Node.ELEMENT_NODE:
                text += SbmlParser._extractTextFromSubtree(node)
            else:
                try:
                    s = node.data.strip()
                    if s:
                        text += s + '\n'
                except AttributeError:
                    pass
        return text


    def _readKineticLaw(self, kineticLawNode):
        """ read subtree under a reaction's 'kineticLaw' node

        Returns True if successful, else False
        """
        isListOfParametersFound = False
        success = True
        reactionId = self._tmpReaction.reactionId

        # kineticLaw has two possible subnodes, math and listOfParameters.
        # We ignore the math node and look for definitions of lower bound,
        # upper bound, flux value, etc. in the listOfParameters.
        for node in kineticLawNode.childNodes:
            nodeName = node.nodeName.lower()

            if nodeName == SbmlDef.LISTOFPARAMETERS.lower():
                if isListOfParametersFound:
                    # If the reaction contains more than one listOfParameters,
                    # only the first occurrence is processed.
                    self._warnDuplicateNode(SbmlDef.LISTOFPARAMETERS,
                                            "reaction", reactionId)
                else:
                    isListOfParametersFound = True
                    success = self._readListOfParameters(node)

            elif (node.nodeType == dom.Node.ELEMENT_NODE and nodeName != ""
                  and nodeName != SbmlDef.MATH.lower()):
                # Report anything that is not math or listOfParameters
                self._debugIgnoringNode(node.nodeName, "reaction", reactionId)

        return success


    def _readListOfParameters(self, listOfParametersNode):
        """ read subtree under the 'listOfParameters' node of a 'kineticLaw'

        This function does not deal with standard SBML, but rather with a
        specific extension for flux balance models.

        Returns True if successful, else False
        """
        isLowerBoundFound = False
        isLowerBoundRead = False
        isUpperBoundFound = False
        isUpperBoundRead = False
        isFluxValFound = False
        isFluxValRead = False
        isObjCoefFound = False
        isObjCoefRead = False
        reactionId = self._tmpReaction.reactionId

        lbDefString = SbmlDef.P_LOWER_BOUND.lower()
        ubDefString = SbmlDef.P_UPPER_BOUND.lower()
        fluxDefString = SbmlDef.P_FLUX_VALUE.lower()
        objCoefDefString = SbmlDef.P_OBJ_COEF.lower()

        for node in listOfParametersNode.childNodes:
            if node.nodeName.lower() == SbmlDef.PARAMETER.lower():
                nodeIdOrig = node.getAttribute(SbmlDef.ID)
                nodeName = node.getAttribute(SbmlDef.NAME).lower()
                nodeId = nodeIdOrig.lower()

                if nodeId == lbDefString or nodeName == lbDefString:
                    if isLowerBoundFound:
                        self._warnDuplicateNode(SbmlDef.P_LOWER_BOUND,
                                                "reaction", reactionId)
                    else:
                        isLowerBoundFound = True
                        if node.getAttribute(SbmlDef.VALUE) != "":
                            try:
                                self._tmpReaction.lb = \
                                    float(node.getAttribute(SbmlDef.VALUE))
                                isLowerBoundRead = True
                            except ValueError:
                                self._errorIllegalFloat(SbmlDef.P_LOWER_BOUND,
                                                       "reaction", reactionId)

                            if self._tmpReaction.lb <= -self.infLevel:
                                self._tmpReaction.lb = -inf
                        else:
                            self._errorMissingValue(SbmlDef.P_LOWER_BOUND,
                                                    "reaction", reactionId)

                elif nodeId == ubDefString or nodeName == ubDefString:
                    if isUpperBoundFound:
                        self._warnDuplicateNode(SbmlDef.P_UPPER_BOUND,
                                                "reaction", reactionId)
                    else:
                        isUpperBoundFound = True
                        if node.getAttribute(SbmlDef.VALUE) != "":
                            try:
                                self._tmpReaction.ub = \
                                    float(node.getAttribute(SbmlDef.VALUE))
                                isUpperBoundRead = True
                            except ValueError:
                                self._errorIllegalFloat(
                                                       SbmlDef.P_UPPER_BOUND,
                                                       "reaction", reactionId)

                            if self._tmpReaction.ub >= self.infLevel:
                                self._tmpReaction.ub = inf
                        else:
                            self._errorMissingValue(SbmlDef.P_UPPER_BOUND,
                                                    "reaction", reactionId)

                elif nodeId == fluxDefString or nodeName == fluxDefString:
                    if isFluxValFound:
                        self._warnDuplicateNode(SbmlDef.P_FLUX_VALUE,
                                                "reaction", reactionId)
                    else:
                        isFluxValFound = True
                        if node.getAttribute(SbmlDef.VALUE) != "":
                            try:
                                self._tmpReaction.fluxVal = \
                                    float(node.getAttribute(SbmlDef.VALUE))
                                isFluxValRead = True
                            except ValueError:
                                self._errorIllegalFloat(SbmlDef.P_FLUX_VALUE,
                                                       "reaction", reactionId)
                        else:
                            self._errorMissingValue(SbmlDef.P_FLUX_VALUE,
                                                    "reaction", reactionId)

                elif nodeId == objCoefDefString or nodeName == objCoefDefString:
                    if isObjCoefFound:
                        self._warnDuplicateNode(SbmlDef.P_OBJ_COEF,
                                                "reaction", reactionId)
                    else:
                        isObjCoefFound = True
                        if node.getAttribute(SbmlDef.VALUE) != "":
                            try:
                                self._tmpReaction.objCoef = \
                                    float(node.getAttribute(SbmlDef.VALUE))
                                isObjCoefRead = True
                            except ValueError:
                                self._errorIllegalFloat(SbmlDef.P_OBJ_COEF,
                                                       "reaction", reactionId)
                        else:
                            self._errorMissingValue(SbmlDef.P_OBJ_COEF,
                                                    "reaction", reactionId)

                elif nodeId == "":
                    msg = ("Warning: Found parameter without attribute \"%s\" "
                           "in %s of reaction '%s'." %
                           (SbmlDef.ID, SbmlDef.LISTOFPARAMETERS, reactionId))
                    self.messages.append((Verbosity.WARNING, msg))

                else:
                    msg = ("Debug: Ignoring parameter %s in %s of reaction '%s'"
                           "." % (nodeIdOrig, SbmlDef.LISTOFPARAMETERS,
                                  reactionId))
                    self.messages.append((Verbosity.DEBUG, msg))

            elif node.nodeType == dom.Node.ELEMENT_NODE and node.nodeName != "":
                # Anything that is not a parameter is ignored
                self._debugIgnoringNode(node.nodeName, "%s of reaction '%s'" %
                                        (SbmlDef.LISTOFPARAMETERS, reactionId))

        if not isLowerBoundFound:
            msg = ("Debug: No %s given for reaction '%s'. Set to default: "
                   "-inf." % (SbmlDef.P_LOWER_BOUND, reactionId))
            self.messages.append((Verbosity.DEBUG, msg))
        if not isUpperBoundFound:
            msg = ("Debug: No %s given for reaction '%s'. Set to default: "
                   "inf." % (SbmlDef.P_UPPER_BOUND, reactionId))
            self.messages.append((Verbosity.DEBUG, msg))

        success = ((not isLowerBoundFound or isLowerBoundRead) and
                   (not isUpperBoundFound or isUpperBoundRead) and
                   (not isFluxValFound or isFluxValRead) and
                   (not isObjCoefFound or isObjCoefRead))

        if (self._tmpReaction.lb > self._tmpReaction.ub):
            msg = ("Warning: %s (%r) > %s (%r) in reaction '%s'. If not fixed, "
                   "this may lead to undefined behavior in simulations." %
                   (SbmlDef.P_LOWER_BOUND, self._tmpReaction.lb,
                    SbmlDef.P_UPPER_BOUND, self._tmpReaction.ub, reactionId))
            self.messages.append((Verbosity.WARNING, msg))

        elif isFluxValRead:
            if (self._tmpReaction.fluxVal < self._tmpReaction.lb or
                self._tmpReaction.fluxVal > self._tmpReaction.ub):
                msg = ("Warning: %s of reaction '%s' (%r) is out of bounds "
                       "(%r, %r)." % (SbmlDef.P_FLUX_VALUE, reactionId,
                       self._tmpReaction.fluxVal, self._tmpReaction.lb,
                       self._tmpReaction.ub))
                self.messages.append((Verbosity.WARNING, msg))

        return success


class CacheReaction(object):
    """ Class for all collected information about a reaction
    """
    def __init__(self):
        self.clear()

    def addReactant(self, reactant):
        self.reactants.append(reactant)

    def addProduct(self, product):
        self.products.append(product)

    def clear(self):
        self.reactionName = ""
        self.reactionId = ""
        self.lb = -inf
        self.ub = inf
        self.reversible = 1
        self.reactants = []
        self.products = []
        self.notes = ""
        self.fluxVal = 0.
        self.objCoef = 0.


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "SBML parser\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-f", "--file", dest="sbmlFile", help="SBML FILE to "
                      "be converted to metano format", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile", help="write "
                      "extracted reactions to FILE", metavar="FILE")
    parser.add_option("-p", "--param", dest="outputParamFile", help="write "
                      "extracted parameters to scenario FILE", metavar="FILE")
    parser.add_option("-q", "--flux", dest="outputFluxFile", help="write fluxes"
                      " (if present) to FILE", metavar="FILE")
    parser.add_option("-i", "--inf-above", dest="infLevel", type="float",
                      help="set LEVEL above which a bound is considered "
                      "infinite (default 1000)", metavar="LEVEL")
    parser.add_option("-n", "--usenames", dest="usenames",
                      help="use metabolite and reaction names instead of IDs",
                      action="store_true")
    parser.add_option("-s", "--suppress-comments", dest="noComments",
                      help="do not write comments to reaction file",
                      action="store_true")
    parser.add_option("-x", "--solver", dest="solver", help="write SOLVER for "
                      "FBA to the scenario file", metavar="SOLVER")
    parser.add_option("-b", "--boundary-regex", dest="boundaryRegex",
                      help="consider all metabolites whose IDs match regular "
                      "expression PATTERN outside the system boundary (in "
                      "addition to metabolites with boundaryCondition=\"true\")"
                      " -- ^xxx matches at beginning, xxx$ at end",
                      metavar="PATTERN")
    parser.add_option("-d", "--debug-info", dest="debugInfo", help="print"
                      " debug information", action="store_true")

    parser.set_defaults(infLevel=1000., debugInfo=False, usenames=False,
                        noComments=False)

    options, _ = parser.parse_args()
    parser.check_required("-f")
    parser.check_required("-p")
    parser.check_required("-o")

    # 2. Set infinity level for bounds

    if options.infLevel <= 0.:
        print("Error: Infinity level must be positive")
        exit()
    sbmlP = SbmlParser(options.infLevel)

    # 3. Read the SBML file

    print("Start parsing..")
    print("SBML file             '%s'" % options.sbmlFile)
    print("Output reaction file: '%s'" % options.outputFile)
    print("Output scenario file: '%s'" % options.outputParamFile)
    if options.boundaryRegex:
        print("Boundary reg-ex:      '%s'" % options.boundaryRegex)
    print("Infinity starts at:   %r\n" % options.infLevel)

    try:
        success = sbmlP.readFile(options.sbmlFile, options.boundaryRegex)
    except ExpatError as strerror:
        print("Error in file %s:" % os.path.basename(options.sbmlFile))
        print(strerror)
        exit()
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.sbmlFile))
        print(strerror)
        exit()

    if options.debugInfo:
        msgs = sbmlP.getMessages(Verbosity.DEBUG)
    else:
        # Show only error, warning, and info messages
        msgs = sbmlP.getMessages()
    if msgs:
        print(('\n'.join(msgs)).encode('ascii', 'replace'))

    # Do not try to build the model if a critical error occurred
    if not success:
        exit()

    # 4. Replace illegal characters (whitespace, :, #) with underscores

    sbmlP.fixNames()

    # 5. Build the model

    model = sbmlP.buildModel(options.usenames)

    if options.usenames:
        # Get unique reaction labels
        reaLabels = dict(list(zip(model.getReactionNames(),
                             model._getTranslatedReactionNames())))
    else:
        reaLabels = None

    # 6. Construct objective function string

    objCoefficients = sbmlP.getObjCoefficients(reaLabels)
    objStr = ""

    if isZeroVector(list(objCoefficients.values())):
        # No coefficients found - look for 'biomass' reaction instead

        # Make a regex pattern for case-insensitive matching of 'biomass'
        pattern = "".join("[%c%c]" % (c, c.lower()) for c in
                          SbmlDef.BIOMASS.upper())
        biomassReactions = [x.name for x in model.getSubmodelByRegex(pattern)]
        nBiomassReactions = len(biomassReactions)
        if nBiomassReactions >= 1:
            objStr = biomassReactions[0]
            if nBiomassReactions > 1:
                print ("Warning: More than one biomass reaction: " +
                       ", ".join(biomassReactions))
        else:
            print ("Warning: No biomass reaction found! (OBJ line must be "
                   "added manually to scenario file.)")
    else:
        # Create objective function string from vector
        for rea in objCoefficients:
            objStr += makeCoefNameTerm(objStr, objCoefficients[rea], rea)

    # 7. Write reaction file

    try:
        model.writeToFile(options.outputFile, not options.noComments)
    except IOError as strerror:
        print ("Unable to write to file %s:" %
               os.path.basename(options.outputFile))
        print(strerror)
        exit()

    # 8. Write scenario file

    try:
        with open(options.outputParamFile,"w") as f:
            preamble = ""
            if objStr:
                preamble += "OBJ max %s\n" % objStr
            if options.solver:
                preamble += "SOLVER %s\n" % options.solver
            f.write(preamble + "\n")
            model.writeBoundsToParamFileHandle(f)
    except IOError as strerror:
        print ("Unable to write to file %s:" %
               os.path.basename(options.outputParamFile))
        print(strerror)
        exit()

    # 9. Write flux file

    if options.outputFluxFile:
        fluxDict = sbmlP.getFluxValues(reaLabels)
        rea_bounds = [(rea.name, (rea.lb, rea.ub)) for rea in model]
        if options.usenames:
            # Replace IDs with names in boundsDict
            boundsDict = {}
            for rea, bounds in rea_bounds:
                boundsDict[reaLabels[rea]] = bounds
        else:
            boundsDict = dict(rea_bounds)

        flux = MetabolicFlux(fluxDict, boundsDict)
        try:
            flux.writeToFile(options.outputFluxFile)
        except IOError as strerror:
            print ("Unable to write to file %s:" %
                   os.path.basename(options.outputFluxFile))
            print(strerror)
            exit()

    print("Finished successfully.")


if __name__ == "__main__":
    main()
