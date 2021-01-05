"""
This module defines constants for structural elements of SBML files.


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

from builtins import object
class SbmlDef(object):
    """ Class for SBML structure elements
        (with extensions for flux balance models)
    """
    # Ubiquitous attributes
    ID = "id"
    NAME = "name"
    NOTES = "notes"
    # "html"/"body" and nested "p"/"html:p" nodes in "notes" are ignored

    # "sbml" nodes and attributes
    SBML = "sbml"
    LEVEL = "level"
    VERSION = "version"
    XMLNS = "xmlns"

    # "model" and child nodes
    MODEL = "model"
    LISTOFUNITDEFINITIONS = "listOfUnitDefinitions"
    LISTOFCOMPARTMENTS = "listOfCompartments"
    LISTOFSPECIES = "listOfSpecies"
    LISTOFREACTIONS = "listOfReactions"

    # listOfUnitDefinitions
    UNITDEFINITION = "unitDefinition"
    LISTOFUNITS = "listOfUnits"
    UNIT = "unit"
    KIND = "kind"
    SCALE = "scale"
    EXPONENT = "exponent"
    MULTIPLIER = "multiplier"

    # listOfCompartments
    COMPARTMENT = "compartment"
    COMPARTMENTOUTSIDE = "outside"

    # listOfSpecies
    SPECIES = "species"
    # "charge" is ignored (and obsolete)
    BOUNDARYCONDITION = "boundaryCondition"
    CONSTANT = "constant"  # can also appear as attribute of speciesReference and compartment
    HASONLYSUBSTANCEUNITS = "hasOnlySubstanceUnits"

    # listOfReactions
    REACTION = "reaction"
    REVERSIBLE = "reversible"
    FAST = "fast"
    LISTOFREACTANTS = "listOfReactants"
    LISTOFPRODUCTS = "listOfProducts"
    KINETICLAW = "kineticLaw"
    MATH = "math"
    CI = "ci"
    CN = "cn"
    LISTOFPARAMETERS = "listOfLocalParameters"

    # listOfReactants / listOfProducts
    SPECIESREFERENCE = "speciesReference"
    STOICHIOMETRY = "stoichiometry"

    # listOfParameters
    PARAMETER = "localParameter"
    VALUE = "value"
    UNITS = "units"
    # Extensions for flux balance models
    P_LOWER_BOUND = "LOWER_BOUND"
    P_UPPER_BOUND = "UPPER_BOUND"
    P_OBJ_COEF = "OBJECTIVE_COEFFICIENT"
    P_FLUX_VALUE = "FLUX_VALUE"

    BIOMASS = "biomass"  # used by SbmlParser to recognize the biomass reaction


def getFirstElementChild(node):
    """ return the first XML child node that is of type Element
    """
    for child in node.childNodes:
        if child.nodeType == node.ELEMENT_NODE:
            return child
