#!/usr/bin/env python
"""
This program exports SBML files from models given in metano's native format. The
exported models conform to SBML Level 2 Version 4, with the Flux Balance
Constraints extension: http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/Flux_Balance_Constraints_(flux).


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

from builtins import str
from builtins import zip
from builtins import range
from builtins import object
from metano.defines import Verbosity, COPYRIGHT_VERSION_STRING
from metano.sbmlstructure import SbmlDef
from metano.fba import OptionParser
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import Reactant, MetabolicModel
from metano.metabolicflux import MetabolicFlux
from numpy import inf
from copy import deepcopy
import xml.dom.minidom as dom
import os
import re


class SbmlExportStruc(object):
    SBML_URL = "http://www.sbml.org/sbml/level3/version1/core"
    XHTML_URL = "http://www.w3.org/1999/xhtml"
    MATH_URL = "http://www.w3.org/1998/Math/MathML"
    XMLNS_HTML = "xmlns:html"
    HTML_P = "html:p"
    ENCODING = "UTF-8"
    DEFAULT_COMPARTMENT = "Cytosol"
    BOUNDARY_COMPARTMENT = "Boundary"
    DEFAULT_BOUNDARY_SUFFIX = "_b"
    THE_FLUX_UNIT = "mmol_per_gDW_per_hr"
    UNITS = ({SbmlDef.KIND: "mole",   SbmlDef.EXPONENT: "1",  SbmlDef.SCALE: "-3", SbmlDef.MULTIPLIER: "1"},
             {SbmlDef.KIND: "gram",   SbmlDef.EXPONENT: "-1",
                 SbmlDef.SCALE: "1",  SbmlDef.MULTIPLIER: "1"},
             {SbmlDef.KIND: "second", SbmlDef.EXPONENT: "-1", SbmlDef.SCALE: "1",  SbmlDef.MULTIPLIER: "0.00027777"})


validIdCharactersPattern = re.compile("[\W_]+")


def removeInvalidIdCharacters(s):

    # remove all special chracters
    s = validIdCharactersPattern.sub("_", s)

    # append '_' to front if ID starts with a number
    s = s if s[0].isalpha() or s[0] == "_" else "_%s" % s

    return s


class SbmlWriter(object):
    """ Class for an SBML writer

    Member variables:

    modelName             -- name of SBML model to be written
    compartments          -- list of compartment IDs
    reactions             -- list of Reaction objects
    reactionLabels        -- list of reaction labels, indexed like reactions
    fluxes                -- flux values (MetabolicFlux object)
    metabolites           -- list of metabolite names
    metaboliteLabels      -- list of metabolite labels, indexed like metabolites
    metaboliteCompartment -- dict { metabolite : compartment }
    isBoundary            -- dict { metabolite : bool } - True for boundary met.
    infinityVal           -- finite value to be used instead of infinite bounds
                             (default: None)
    messages              -- error, warning, info, and debug messages generated
                             while exporting the model
    """

    _encodeBool = {False: "false", True: "true"}

    def __init__(self, model, modelName, maxmin, objStr="",
                 fluxes=MetabolicFlux(), stripComments=False, infinityVal=None):
        """ initialize class SMBLParser

        Keyword arguments:

        model              -- MetabolicModel to be exported to SBML
        modelName          -- name of SBML model to be written
        maxmin             -- True for maximization, False for minimization
        objStr             -- objective function definition (given as string)
        fluxes             -- flux values (MetabolicFlux object)
        infinityVal        -- finite value to be used instead of infinite bounds
                              (default: None)
        """
        self.clear()
        self.modelName = modelName

        useLabelsBak = model.useReactionLabels, model.useMetaboliteLabels
        model.useReactionLabels = True
        model.useMetaboliteLabels = True
        self.reactions = deepcopy(model.reactions)
        self.reactionLabels = deepcopy(model._getTranslatedReactionNames())
        self.reactionIndex = deepcopy(model.reactionDict)
        self.metabolites = deepcopy(model.metabolites)
        self.metaboliteLabels = deepcopy(model._getTranslatedMetaboliteNames())
        self.metaboliteIndex = deepcopy(model.metaboliteDict)
        model.useReactionLabels, model.useMetaboliteLabels = useLabelsBak

        self.isBoundary = dict((met, False) for met in self.metabolites)

        # Set coefficient vector of objective function
        try:
            # maxmin is inverted because in our SBML files (unlike in LP),
            # maximization is the default (i.e. factor 1 instead of -1)
            self.objCoef = ParamParser.convertObjFuncToLinVec(objStr,
                                                              model.reactionDict, maxmin=not maxmin)
        except ValueError:
            msg = ("Info: Objective function (%s) is not a linear expression. "
                   "%s will not be written." % (objStr, SbmlDef.P_OBJ_COEF))
            self.messages.append((Verbosity.INFO, msg))
        except KeyError as strerror:
            msg = ("Warning: An error occurred in evaluating the objective "
                   "function:\n %s." % strerror)
            self.messages.append((Verbosity.WARNING, msg))

        self.fluxes = fluxes
        self.stripComments = stripComments

        if infinityVal and infinityVal < 0.:
            infinityVal = -1.*infinityVal
        self.infinityVal = infinityVal

    def clear(self):
        """ clear all internal variables, warnings and info messages
        """
        self.modelName = ""
        self.compartments = []
        self.reactions = []
        self.reactionLabels = []
        self.fluxes = MetabolicFlux()
        self.objCoef = []
        self.metabolites = []
        self.metaboliteLabels = []
        self.metaboliteCompartment = {}
        self.isBoundary = {}
        self.messages = []

    def getMessages(self, level=Verbosity.INFO):
        """ return a list of all messages at or above the given level of
            severity

        Levels are defined in defines.py.
        """
        return [x[1] for x in self.messages if x[0] <= level]

    def makeCompartments(self, compartmentRegex="",
                         defaultCompartment=SbmlExportStruc.DEFAULT_COMPARTMENT):
        """ construct compartment IDs and metabolite->compartment mapping from
            self.metabolites

        Keyword arguments:

        compartmentRegex   -- regular expression describing what part (if any)
                              of a metabolite name identifies the compartment -
                              must contain exactly one group, which is used to
                              generate compartment IDs
        defaultCompartment -- ID of compartment of metabolites not matching
                              compartmentRegex

        Modified member variables:

        compartments, metaboliteCompartment

        Returns True if successful, else False
        """
        if not compartmentRegex:
            # Assign all metabolites to the same (default) compartment
            self.compartments = [defaultCompartment]
            self.metaboliteCompartment = dict((met, defaultCompartment)
                                              for met in self.metabolites)
            msg = ("Info: No compartment reg-ex given. Assigned all metabolites"
                   " to default compartment (%s)." % defaultCompartment)
            self.messages.append((Verbosity.INFO, msg))
            return True

        try:
            pattern = re.compile(compartmentRegex)
        except re.error as strerror:
            msg = ("Error in regular expression: %s." %
                   strerror)
            self.messages.append((Verbosity.ERROR, msg))
            return False

        setOfCompartments = set()

        for met in self.metabolites:
            reHit = pattern.search(met)
            if not reHit:
                # Assign non-matching metabolites to default compartment
                compartment = defaultCompartment
                msg = ("Debug: Metabolite '%s' does not match compartment "
                       "reg-ex. Assigned to default compartment (%s)." %
                       (met, defaultCompartment))
                self.messages.append((Verbosity.DEBUG, msg))
            else:
                try:
                    compartment = reHit.group(1)
                except IndexError:
                    msg = ("Error: Regular expression has no groups. Use "
                           "parentheses to define one.\n"
                           "       Example: \"_(ex|c|p)$\".")
                    self.messages.append((Verbosity.ERROR, msg))
                    return False

            compartment = removeInvalidIdCharacters(compartment)
            self.metaboliteCompartment[met] = compartment

            if compartment not in setOfCompartments:
                self.compartments.append(compartment)
                setOfCompartments.add(compartment)
                msg = "Info: Created compartment '%s'." % compartment
                self.messages.append((Verbosity.INFO, msg))

        return True

    def addBoundaryMetabolites(self,
                               suffix=SbmlExportStruc.DEFAULT_BOUNDARY_SUFFIX):
        """ add boundary metabolites to all reactions that have an empty left-
            or right-hand side

        Boundary metabolites are marked with the given suffix and assigned their
        own compartment ("Boundary"). No non-boundary metabolite may have the
        designated boundary suffix.

        Returns True if successful, else False
        """
        for metabolite in self.metabolites:
            if metabolite.endswith(suffix):
                msg = ("Error: Non-boundary metabolite '%s' found with "
                       "designated boundary suffix '%s'. Please choose a unique"
                       " suffix." % (metabolite, suffix))
                self.messages.append((Verbosity.ERROR, msg))
                return False

        hasBoundaryCompartment = False
        incompleteReactions = [r for r in self.reactions if not r.getEducts() or
                               not r.getProducts()]

        for reaction in incompleteReactions:
            educts = reaction.getEducts()
            products = reaction.getProducts()

            if not educts and not products:
                # Empty reaction, illegal in SBML (also in metano)
                i = self.reactionIndex[reaction.name]
                del self.reactions[i]
                del self.reactionLabels[i]
                # Rebuild index
                self.reactionIndex = dict(list(zip([r.name for r in self.reactions],
                                                   list(range(len(self.reactions))))))
                msg = ("Warning: Reaction '%s' has neither reactants "
                       "nor products. Deleted." % reaction.name)
                self.messages.append((Verbosity.WARNING, msg))

            elif products == []:
                if not hasBoundaryCompartment:
                    hasBoundaryCompartment = self._addBoundaryCompartment()

                self._addBoundaryMetabolitesToSide(reaction, "right", educts,
                                                   suffix)
            elif educts == []:
                if not hasBoundaryCompartment:
                    hasBoundaryCompartment = self._addBoundaryCompartment()

                self._addBoundaryMetabolitesToSide(reaction, "left", products,
                                                   suffix)
        return True

    def _addBoundaryCompartment(self):
        """ add "Boundary" compartment to list of compartments (if not already
            present)

            Returns True
        """
        if SbmlExportStruc.BOUNDARY_COMPARTMENT not in self.compartments:
            self.compartments.append(SbmlExportStruc.BOUNDARY_COMPARTMENT)
        return True

    def _addBoundaryMetabolitesToSide(self, reaction, side, reactants, suffix):
        """ add boundary metabolites to the specified side of the given reaction

        Keyword arguments:

        reaction  -- Reaction object to be modified
        side      -- "left" or "right"
        reactants -- list of Reactant objects for the other (non-boundary) side
        suffix    -- suffix for marking the newly created metabolites
        """
        for reactant in reactants:
            boundaryName = reactant.name + suffix
            boundaryReactant = Reactant(boundaryName, -1.*reactant.coef)

            if boundaryName not in self.metaboliteCompartment:
                self.metabolites.append(boundaryName)
                self.metaboliteCompartment[boundaryName] = \
                    SbmlExportStruc.BOUNDARY_COMPARTMENT
                self.isBoundary[boundaryName] = True
                orig_index = self.metaboliteIndex[reactant.name]
                boundaryLabel = (self.metaboliteLabels[orig_index] + suffix)
                self.metaboliteLabels.append(boundaryLabel)
                self.metaboliteIndex[boundaryName] = len(self.metaboliteIndex)

            reaction.addReactant(boundaryReactant)
            msg = ("Debug: Added boundary species '%s' with "
                   "coefficient %r to %s-hand side of reaction '%s'."
                   % (boundaryName, abs(reactant.coef), side, reaction.name))
            self.messages.append((Verbosity.DEBUG, msg))

    def writeFile(self, sbmlFileName):
        """ write prepared structures to file designated by sbmlFileName
        """
        sbmlTree = dom.Document()

        sbml = dom.Element(SbmlDef.SBML)
        sbml.setAttribute(SbmlDef.LEVEL, "3")
        sbml.setAttribute(SbmlDef.VERSION, "1")
        sbml.setAttribute(SbmlDef.XMLNS, SbmlExportStruc.SBML_URL)
        sbml.setAttribute(SbmlExportStruc.XMLNS_HTML,
                          SbmlExportStruc.XHTML_URL)

        model = dom.Element(SbmlDef.MODEL)
        model.setAttribute(
            SbmlDef.ID, removeInvalidIdCharacters(self.modelName))
        model.setAttribute(SbmlDef.NAME, self.modelName)

        # Write unit definition
        listOfUnitDefinitions = dom.Element(SbmlDef.LISTOFUNITDEFINITIONS)
        unitDefinition = dom.Element(SbmlDef.UNITDEFINITION)
        unitDefinition.setAttribute(
            SbmlDef.ID, removeInvalidIdCharacters(SbmlExportStruc.THE_FLUX_UNIT))
        listOfUnits = dom.Element(SbmlDef.LISTOFUNITS)
        for u in SbmlExportStruc.UNITS:
            unit = dom.Element(SbmlDef.UNIT)
            for attr in u:
                unit.setAttribute(attr, u[attr])
            listOfUnits.appendChild(unit)
        unitDefinition.appendChild(listOfUnits)
        listOfUnitDefinitions.appendChild(unitDefinition)
        model.appendChild(listOfUnitDefinitions)

        # Write compartments
        listOfCompartments = dom.Element(SbmlDef.LISTOFCOMPARTMENTS)
        for compartment in self.compartments:
            element = dom.Element(SbmlDef.COMPARTMENT)
            element.setAttribute(
                SbmlDef.ID, removeInvalidIdCharacters(compartment))
            element.setAttribute(SbmlDef.CONSTANT, "true")
            listOfCompartments.appendChild(element)
        model.appendChild(listOfCompartments)

        # Create list of species
        listOfSpecies = dom.Element(SbmlDef.LISTOFSPECIES)
        for i in range(len(self.metabolites)):
            metabolite = self.metabolites[i]
            metaboliteLabel = self.metaboliteLabels[i]
            metaboliteCompartment = self.metaboliteCompartment[metabolite]
            boundaryCondition = self._encodeBool[self.isBoundary[metabolite]]
            species = dom.Element(SbmlDef.SPECIES)
            species.setAttribute(
                SbmlDef.ID, removeInvalidIdCharacters(metabolite))
            species.setAttribute(SbmlDef.NAME, metaboliteLabel)
            species.setAttribute(SbmlDef.COMPARTMENT, metaboliteCompartment)
            species.setAttribute(SbmlDef.BOUNDARYCONDITION, boundaryCondition)
            species.setAttribute(SbmlDef.HASONLYSUBSTANCEUNITS, "true")
            species.setAttribute(SbmlDef.CONSTANT, "false")
            listOfSpecies.appendChild(species)

            msg = "Debug: Created species node '%s' " % metabolite
            if metaboliteLabel and metaboliteLabel != metabolite:
                msg += "(%s) " % metaboliteLabel
            msg += ("with %s='%s' and %s='%s'." % (SbmlDef.COMPARTMENT,
                                                   metaboliteCompartment, SbmlDef.BOUNDARYCONDITION,
                                                   boundaryCondition))
            self.messages.append((Verbosity.DEBUG, msg))

        model.appendChild(listOfSpecies)

        # Create list of reactions
        listOfReactions = self._createListOfReactions()
        model.appendChild(listOfReactions)

        sbml.appendChild(model)
        sbmlTree.appendChild(sbml)

        # Write tree to file
        with open(sbmlFileName, 'w') as f:
            sbmlTree.writexml(f, "", '\t', '\n', SbmlExportStruc.ENCODING)

    def _createListOfReactions(self):
        listOfReactions = dom.Element(SbmlDef.LISTOFREACTIONS)
        fluxVec = self.fluxes.getVecOrderedByDict(self.reactionIndex)
        hasFluxes = len([f for f in fluxVec if f is not None]) != 0
        if self.fluxes and not hasFluxes:
            # Generate warning if model and flux distribution do not share any
            # reaction (else: generate a warning for each reaction without flux)
            msg = ("Warning: The given flux distribution does not fit the model"
                   " (reaction names are different).")
            self.messages.append((Verbosity.WARNING, msg))

        for i in range(len(self.reactions)):
            # Set required attributes for reaction
            rea_item = self.reactions[i]
            reactionName = rea_item.name
            reactionLabel = self.reactionLabels[i]
            reversible = self._encodeBool[rea_item.direction == 0]
            reaction = dom.Element(SbmlDef.REACTION)
            reaction.setAttribute(
                SbmlDef.ID, removeInvalidIdCharacters(reactionName))
            reaction.setAttribute(SbmlDef.NAME, reactionLabel)
            reaction.setAttribute(SbmlDef.REVERSIBLE, reversible)
            reaction.setAttribute(SbmlDef.FAST, "false")

            msg = "Debug: Created reaction node '%s' " % reactionName
            if reactionLabel and reactionLabel != reactionName:
                msg += "(%s) " % reactionLabel
            msg += "with %s='%s'." % (SbmlDef.REVERSIBLE, reversible)
            self.messages.append((Verbosity.DEBUG, msg))

            if rea_item.direction == -1:
                # Irreversible reaction proceeding from right to left must be
                # reversed - in SBML, irreversible implies left to right
                educts, products = rea_item.getProducts(), rea_item.getEducts()
                tmp_ub, tmp_lb = -rea_item._lb, -rea_item._ub
                tmp_flux = -fluxVec[i] if fluxVec[i] is not None else None
                msg = ("Debug: Replacing irreversible right-to-left reaction "
                       "'%s' with left-to-right reaction." % reactionName)
                self.messages.append((Verbosity.DEBUG, msg))
            else:
                educts, products = rea_item.getEducts(), rea_item.getProducts()
                tmp_lb, tmp_ub = rea_item.getBounds()
                tmp_flux = fluxVec[i]

            # Create notes
            if not self.stripComments:
                comments = str(rea_item.comment).splitlines()
                if len(comments) > 0:
                    notes = dom.Element(SbmlDef.NOTES)
                    for line in comments:
                        paragraph = dom.Element(SbmlExportStruc.HTML_P)
                        t = dom.Text()
                        t.data = line.strip()
                        paragraph.appendChild(t)
                        notes.appendChild(paragraph)
                    reaction.appendChild(notes)

            # Create list of reactants and list of products
            if educts:
                listOfReactants = self._createListOfReactantsProducts(
                    SbmlDef.LISTOFREACTANTS, educts)
                reaction.appendChild(listOfReactants)
            if products:
                listOfProducts = self._createListOfReactantsProducts(
                    SbmlDef.LISTOFPRODUCTS, products)
                reaction.appendChild(listOfProducts)

            # Create kineticLaw for reaction
            kineticLaw = dom.Element(SbmlDef.KINETICLAW)

            math = dom.Element(SbmlDef.MATH)
            math.setAttribute(SbmlDef.XMLNS, SbmlExportStruc.MATH_URL)
            # ===================================================================
            # ci = dom.Element(SbmlDef.CI)
            # t = dom.Text()
            # t.data = "%s" % SbmlDef.P_FLUX_VALUE
            # ci.appendChild(t)
            # math.appendChild(ci)
            # ===================================================================
            cn = dom.Element(SbmlDef.CN)
            t = dom.Text()

            # Write listOfParameters with lower bound, upper bound, objective
            # coefficient, and flux value
            listOfParameters = dom.Element(SbmlDef.LISTOFPARAMETERS)
            paramMsg = "Debug: Parameters: "

            if self.infinityVal or tmp_lb != -inf:

                if tmp_lb != -inf:
                    lb = tmp_lb
                else:
                    lb = -1*self.infinityVal

                parameter = dom.Element(SbmlDef.PARAMETER)
                parameter.setAttribute(
                    SbmlDef.ID, removeInvalidIdCharacters(SbmlDef.P_LOWER_BOUND))
                # parameter.setAttribute(SbmlDef.CONSTANT, "true") # Change in SBML V3
                value = repr(lb)
                parameter.setAttribute(SbmlDef.VALUE, value)
                parameter.setAttribute(SbmlDef.UNITS,
                                       SbmlExportStruc.THE_FLUX_UNIT)
                listOfParameters.appendChild(parameter)
                paramMsg += "%s: %r" % (SbmlDef.P_LOWER_BOUND, lb)
            else:
                paramMsg += "no lower bound"

            if self.infinityVal or tmp_ub != inf:

                if tmp_ub != inf:
                    ub = tmp_ub
                else:
                    ub = self.infinityVal

                parameter = dom.Element(SbmlDef.PARAMETER)
                parameter.setAttribute(
                    SbmlDef.ID, removeInvalidIdCharacters(SbmlDef.P_UPPER_BOUND))
                # parameter.setAttribute(SbmlDef.CONSTANT, "true") # Change in SBML V3
                value = repr(ub)
                parameter.setAttribute(SbmlDef.VALUE, value)
                parameter.setAttribute(SbmlDef.UNITS,
                                       SbmlExportStruc.THE_FLUX_UNIT)
                listOfParameters.appendChild(parameter)
                paramMsg += ", %s: %r" % (SbmlDef.P_UPPER_BOUND, ub)
            else:
                paramMsg += ", no upper bound"

            if self.objCoef:
                coef = self.objCoef[i]
                parameter = dom.Element(SbmlDef.PARAMETER)
                parameter.setAttribute(
                    SbmlDef.ID, removeInvalidIdCharacters(SbmlDef.P_OBJ_COEF))
                # parameter.setAttribute(SbmlDef.CONSTANT, "true") # Change in SBML V3
                value = repr(coef)
                parameter.setAttribute(SbmlDef.VALUE, value)
                listOfParameters.appendChild(parameter)
                paramMsg += ", %s: %r" % (SbmlDef.P_OBJ_COEF, coef)

            if tmp_flux is not None:
                parameter = dom.Element(SbmlDef.PARAMETER)
                parameter.setAttribute(
                    SbmlDef.ID, removeInvalidIdCharacters(SbmlDef.P_FLUX_VALUE))
                # parameter.setAttribute(SbmlDef.CONSTANT, "true") # Change in SBML V3
                value = repr(tmp_flux)
                parameter.setAttribute(SbmlDef.VALUE, value)
                parameter.setAttribute(SbmlDef.UNITS,
                                       SbmlExportStruc.THE_FLUX_UNIT)
                listOfParameters.appendChild(parameter)
                paramMsg += ", %s: %r" % (SbmlDef.P_FLUX_VALUE, tmp_flux)
            elif hasFluxes:
                msg = ("Warning: Reaction '%s' has no flux value." %
                       reactionName)
                self.messages.append((Verbosity.WARNING, msg))

            t.data = "%s" % value
            cn.appendChild(t)
            math.appendChild(cn)
            kineticLaw.appendChild(math)

            kineticLaw.appendChild(listOfParameters)
            reaction.appendChild(kineticLaw)
            self.messages.append((Verbosity.DEBUG, paramMsg))

            # Store reaction in the reaction list
            listOfReactions.appendChild(reaction)

        return listOfReactions

    def _createListOfReactantsProducts(self, typ, reactants):
        """ create listOfReactants or listOfProducts node

        Keyword arguments:

        typ       -- SbmlDef.LISTOFREACTANTS or SbmlDef.LISTOFPRODUCTS
        reactants -- list of Reactant objects

        Returns SBML node with the given type
        """
        listOfSpeciesRef = dom.Element(typ)
        for r in reactants:
            reactant = dom.Element(SbmlDef.SPECIESREFERENCE)
            reactant.setAttribute(
                SbmlDef.SPECIES, removeInvalidIdCharacters(r.name))
            reactant.setAttribute(SbmlDef.CONSTANT, "false")
            reactant.setAttribute(SbmlDef.STOICHIOMETRY, repr(abs(r.coef)))
            listOfSpeciesRef.appendChild(reactant)

        return listOfSpeciesRef


def main():
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "SBML writer\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)

    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="write an SBML file for the network given "
                           "by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="use the given scenario FILE for SBML export",
                      metavar="FILE")
    parser.add_option("-q", "--flux", dest="fluxFile", help="write fluxes read "
                      "from FILE to SBML file (optional)", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write output to SBML FILE",
                      metavar="FILE")
    parser.add_option("-c", "--compartment-regex", dest="compartmentRegex",
                      help="regular expression PATTERN in metabolite names that"
                      " identifies compartment information (e.g. "
                      "\"\\[(\\w+)\\]$\" or \"_(ex|c|p)$\") -- must contain "
                      "exactly one group (in parentheses)", metavar="PATTERN")
    parser.add_option("-d", "--default-compartment", dest="defaultCompartment",
                      help="NAME of default compartment (assigned to all "
                      "metabolites not matching the compartment reg-ex, "
                      "default: Cytosol)", metavar="NAME")
    parser.add_option("-s", "--suppress-comments", dest="noComments",
                      help="do not transfer comments to SBML file",
                      action="store_true")
    parser.add_option("-i", "--infinity", dest="infinityVal", type="float",
                      help="use finite VALUE instead of infinity (defaults to "
                      "None)", metavar="VALUE")
    parser.add_option("-e", "--debug-info", dest="debugInfo", help="print"
                      " debug information", action="store_true")
    parser.add_option("-b", "--boundary-suffix", dest="boundarySuffix",
                      help="add boundary metabolites, marked with the given "
                      "SUFFIX to reactions lacking either products or "
                      "reactants", metavar="SUFFIX")

    parser.set_defaults(defaultCompartment=SbmlExportStruc.DEFAULT_COMPARTMENT,
                        noComments=False, debugInfo=False)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-o")

    print("Start creating SBML file '%s'.." % options.outputFile)
    print("Reaction file:           '%s'" % options.reactionFile)
    print("Scenario file:           '%s'" % options.paramFile)
    if options.fluxFile:
        print("Flux file:               '%s'" % options.fluxFile)
    print("Compartment reg-ex:      '%s'" % options.compartmentRegex)
    print("Default compartment:     '%s'" % options.defaultCompartment)
    print("Infinity value:          %r" % options.infinityVal)
    if options.boundarySuffix:
        print("Boundary suffix:         '%s'\n" % options.boundarySuffix)

    # 2. Parse reaction file

    rparser = ReactionParser()
    model = MetabolicModel()
    try:
        model.addReactionsFromFile(options.reactionFile, rparser)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.reactionFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in reaction file %s:" %
               os.path.basename(options.reactionFile))
        print(strerror)
        exit()

    # 3. Parse scenario file

    model_messages = []
    pparser = ParamParser()
    try:
        # Read objective function and bounds from scenario file
        maxmin, objStr, _, _, lb, ub = pparser.parse(options.paramFile)
        # Set flux bounds in model
        model.setFiniteBounds(lb, ub, True, model_messages)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(options.paramFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in scenario file %s:" %
               os.path.basename(options.paramFile))
        print(strerror)
        exit()
    except ValueError as strerror:
        print(strerror)
        exit()

    # Show warning and info messages of parsers
    msgs = (rparser.getMessages() + pparser.getMessages() +
            [x[1] for x in model_messages])
    if msgs:
        print('\n'.join(msgs))

    # 4. Read flux distribution from file

    fluxes = MetabolicFlux()
    if options.fluxFile:
        try:
            fluxes.readFromFile(options.fluxFile)
        except IOError as strerror:
            print ("An error occurred while trying to read file %s:" %
                   os.path.basename(options.fluxFile))
            print(strerror)
            print("Fluxes will not be written.")
            fluxes = MetabolicFlux()
        except SyntaxError as strerror:
            print ("An error occurred parsing file %s:" %
                   os.path.basename(options.fluxFile))
            print(strerror)
            print("Fluxes will not be written.")
            fluxes = MetabolicFlux()

    # 5. Write SBML file

    # Name of model is name of output file without file extension
    modelName = os.path.basename(options.outputFile)
    dotPos = modelName.rfind(".")
    if dotPos > 0:
        modelName = modelName[:dotPos]

    sbmlWriter = SbmlWriter(model, modelName, maxmin, objStr, fluxes,
                            options.noComments, options.infinityVal)
    success = sbmlWriter.makeCompartments(options.compartmentRegex,
                                          options.defaultCompartment)
    if success:
        if options.boundarySuffix:
            success = sbmlWriter.addBoundaryMetabolites(options.boundarySuffix)
        if success:
            sbmlWriter.writeFile(options.outputFile)

    if options.debugInfo:
        msgs = sbmlWriter.getMessages(Verbosity.DEBUG)
    else:
        # Show only error, warning, and info messages
        msgs = sbmlWriter.getMessages()
    if msgs:
        print('\n'+'\n'.join(msgs))

    print("Finished.")


if __name__ == "__main__":
    main()
