"""
This module defines the class ReactionParser, a parser for reaction files.


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

from builtins import map
from builtins import object
from metano.defines import Verbosity, ReaFileStruc


def _splitAtPlus(s):
    """ split the string at plus signs that are preceded or followed by
        whitespace
    """
    tokens = list(map(str.strip, s.split(' +')))
    tokensList = [x.split('\t+') for x in tokens]
    tokens = []
    for t in tokensList:
        tokens += list(map(str.strip, t))
    tokensList = [x.split('+ ') for x in tokens]
    tokens = []
    for t in tokensList:
        tokens += list(map(str.strip, t))
    tokensList = [x.split('+\t') for x in tokens]
    tokens = []
    for t in tokensList:
        tokens += list(map(str.strip, t))
    return tokens


class ReactionParser(object):
    """Parser for reaction files.
    It reads lines of the format

    <reaction> : <compound> [+ <compound>...] -->|<--|<=>
                 <compound> [+ <compound>...]

    and builds the stoichiometric matrix of the corresponding metabolic network,
    along with dictionaries of reactions and compounds, which relate names to
    the corresponding rows or columns in the matrix, and a list of reversibility
    flags for all reactions.
    """

    format_hint = ("Format is\n<reaction> %s <compound> [+ <compound>...] "
                   "%s|%s|%s <compound> [+ <compound>...]" %
                   (ReaFileStruc.delimiter, ReaFileStruc.arrowIrr,
                    ReaFileStruc.arrowFlip, ReaFileStruc.arrowRev))

    def __init__(self, skipCommentLines=False):
        """ initialize the ReactionParser
        """
        self.clear()
        self.skipCommentLines = skipCommentLines

    def clear(self):
        """ clear all internal variables, warnings and info messages
        """
        # List of warnings and info messages generated during parsing,
        # elements are pairs (level, msg)
        self.messages = []
        # If true, don't assign rows starting with a hash sign (#) to the
        # reaction below the comment
        self.skipCommentLines = False

    def getMessages(self, level=Verbosity.INFO):
        """ return a list of all messages at or above the given level of
            severity
            Levels are defined in defines.py.
        """
        return [x[1] for x in self.messages if x[0] <= level]

    def parse(self, filename):
        """ parse the reaction file identified by filename
            -- wrapper for parseByHandle()
        """
        with open(filename) as f:
            return self.parseByHandle(f)

    def parseByHandle(self, f):
        """ parse the reaction file given as a file object

        Returns:
        list of (name, educts, arrow, products, comment) tuples with

        name     -- reaction name
        educts   -- list of (coef, name) tuples for left-hand-side metabolites
        arrow    -- reaction arrow (<=>, -->, or <--)
        products -- list of (coef, name) tuples for right-hand-side metabolites
        comment  -- an optional comment for the reaction
        """
        line_no = 0
        result = []
        comment = ""
        metabSet, reactionSet = set(), set()
        for line in f:
            line_no += 1
            if line.startswith(ReaFileStruc.commentSign):
                continue

            # Everything after comment sign ('#') is comment
            comment_pos = line.find(ReaFileStruc.commentSign)
            if comment_pos > 0 or (comment_pos == 0 and
                                   not self.skipCommentLines):
                if comment:
                    # Append to previous line of comment
                    comment += '\n'
                comment += line[comment_pos+1:].strip()
            if comment_pos >= 0:
                line = line[:comment_pos].strip()

            if line == "" or line.isspace():
                if comment and comment_pos < 0:
                    comment += '\n'
                continue  # Skip blank lines
            else:
                # Current line contains a reaction -> finish comment
                reaction_comment = comment
                comment = ""

            # Split line at delimiter
            delimiter_pos = line.find(ReaFileStruc.delimiter)
            if delimiter_pos < 0:
                msg = self.error(line_no)
                raise SyntaxError(msg)

            rea_name = line[:delimiter_pos].strip()
            if rea_name == "":  # Reaction name must not be empty
                msg = self.error(line_no)
                raise SyntaxError(msg)

            if rea_name in reactionSet:
                msg = self.error(line_no, "Duplicate definition of reaction " +
                                 rea_name, False)
                raise SyntaxError(msg)
            if len(rea_name.split(None, 1)) > 1:
                msg = ("Warning: Reaction name '%s' contains space/tab "
                       "(reaction file, line %u)" % (rea_name, line_no))
                self.messages.append((Verbosity.WARNING, msg))
            reactionSet.add(rea_name)

            rest = line[delimiter_pos+len(ReaFileStruc.delimiter):]
            if rest.find(ReaFileStruc.delimiter) >= 0:
                msg = self.error(line_no, "Extra delimiter ('%s')" %
                                 ReaFileStruc.delimiter)
                raise SyntaxError(msg)

            # Split rest at arrow
            for arrow in ReaFileStruc.arrow_prop:
                arrow_pos = rest.find(arrow)
                if arrow_pos >= 0:
                    arrowprop = ReaFileStruc.arrow_prop[arrow]
                    break

            if arrow_pos < 0:
                msg = self.error(line_no)
                raise SyntaxError(msg)

            # Get reactants and products
            educts_str = rest[:arrow_pos]
            products_str = rest[arrow_pos+arrowprop.length:]

            line_minus_arrow = educts_str + '\n' + products_str
            # ('\n' is never part of arrow string)
            for arrow_ in ReaFileStruc.arrow_prop:
                if line_minus_arrow.find(arrow_) >= 0:
                    msg = self.error(line_no, "Extra arrow")
                    raise SyntaxError(msg)

            educts = _splitAtPlus(educts_str)
            products = _splitAtPlus(products_str)

            # Get reactants as (coefficient, name) pairs
            educts = self.process_metabolites(educts, line_no)
            products = self.process_metabolites(products, line_no)

            # Check for duplicate metabolites in reaction
            metabNames = set()
            for name in [x[1] for x in educts+products]:
                if name in metabNames:
                    msg = ("Warning: Duplicate metabolite '%s' in line %u in "
                           "reaction file." % (name, line_no))
                    self.messages.append((Verbosity.WARNING, msg))
                else:
                    metabNames.add(name)
                    metabSet.add(name)

            result.append((rea_name, educts, arrow,
                           products, reaction_comment))

        msg = ("Info: The metabolic network has %u reactions and %u metabolites"
               ".") % (len(reactionSet), len(metabSet))
        self.messages.append((Verbosity.INFO, msg))
        return result

    def process_metabolites(self, metab_strings, line_no):
        """ process a list of metabolites in a reaction (reactants or products)

        Keyword arguments:

        metab_strings -- list of metabolites given as strings in the form
                             "<coefficient> <metabolite>"
        line_no       -- line number in file

        Returns: list of (coefficient, name) pairs
        """
        result = []
        for metab_string in metab_strings:
            if metab_string == "":  # Skip empty strings
                continue

            entry = metab_string.split(None, 2)
            if len(entry) > 2:
                msg = self.error(line_no, "Metabolite names must not contain "
                                 "spaces/tabs.\n Illegal name encountered")
                raise SyntaxError(msg)
            elif len(entry) == 1:
                name = entry[0]
                coefficient = 1.
            else:
                # Split entry into coefficient and name
                try:
                    coefficient = float(entry[0])
                except ValueError:
                    msg = self.error(line_no, "Metabolite names must not "
                                     "contain spaces/tabs.\n Illegal name "
                                     "encountered", False)
                    raise SyntaxError(msg)
                name = entry[1]

            if coefficient <= 0.:
                msg = self.error(line_no, "Negative or zero coefficient "
                                 "encountered", False)
                raise SyntaxError(msg)

            result.append((coefficient, name))
        return result

    def error(self, line_no=0, msg="Malformed reaction", hint=True):
        """ generate an error message

        Keyword arguments:

        line_no     -- the line number, 0 if not applicable
        msg         -- error message
        hint        -- show hint on expected format of reactions (True/False)
        """
        s = msg
        if line_no > 0:
            s += " in line %u" % line_no
        if hint:
            s += "\n"
            s += self.format_hint
        return s

    def fixfile(self, filename, outfilename):
        """ fix reaction file by replacing spaces in reaction and metabolite
            names with underscores and converting metabolite names to lowercase

        Keyword arguments:

        filename    -- input file
        outfilename -- name of output file to be written
        """
        with open(filename) as f:
            with open(outfilename, 'w') as outf:
                self.fixfileByHandles(f, outf)

    def fixfileByHandles(self, f, outf):
        """ fix reaction file by replacing spaces in reaction and metabolite
            names with underscores and converting metabolite names to lowercase

        Keyword arguments:

        f    -- file handle of input file (must be open for reading)
        outf -- file handle of output file (must be open for writing)
        """
        line_no = 0
        reaNamesFixed = 0
        metabDict = {}  # dictionary {name : new name}
        for line in f:
            line_no += 1

            # Everything after comment sign ('#') is comment
            comment_pos = line.find(ReaFileStruc.commentSign)
            if comment_pos >= 0:
                # Get comment with comment sign
                comment = line[comment_pos:].rstrip()
            else:
                comment = ""
            line = line[:comment_pos]
            if line == "" or line.isspace():
                outf.write(comment+'\n')
                continue  # Skip blank lines

            # Split line at delimiter
            delimiter_pos = line.find(ReaFileStruc.delimiter)
            if delimiter_pos < 0:
                msg = self.error(line_no)
                raise SyntaxError(msg)

            # Fix reaction name
            rea_name = line[:delimiter_pos].strip()
            if rea_name == "":  # Reaction name must not be empty
                msg = self.error(line_no)
                raise SyntaxError(msg)

            rea_name_split = rea_name.split()
            if len(rea_name_split) > 1:
                reaNamesFixed += 1
                rea_name_new = '_'.join(rea_name_split)
                prefix = (rea_name_new + " %s " %
                          ReaFileStruc.delimiter)
                msg = "Reaction '%s' --> '%s'" % (rea_name, rea_name_new)
                self.messages.append((Verbosity.DEBUG, msg))
            else:
                prefix = line[:delimiter_pos+len(ReaFileStruc.delimiter)] + ' '

            rest = line[delimiter_pos+len(ReaFileStruc.delimiter):]
            if rest.find(ReaFileStruc.delimiter) >= 0:
                msg = self.error(line_no, "Extra delimiter ('%s')" %
                                 ReaFileStruc.delimiter)
                raise SyntaxError(msg)

            # Split rest at arrow
            for arrow in ReaFileStruc.arrow_prop:
                arrow_pos = rest.find(arrow)
                if arrow_pos >= 0:
                    arrowprop = ReaFileStruc.arrow_prop[arrow]
                    arrowStr = ' ' + arrow + ' '
                    break

            if arrow_pos < 0:
                msg = self.error(line_no)
                raise SyntaxError(msg)

            # Get reactants and products
            educts_str = rest[:arrow_pos]
            products_str = rest[arrow_pos+arrowprop.length:]

            line_minus_arrow = educts_str + '\n' + products_str
            # ('\n' is never part of arrow string)
            for arrow_ in ReaFileStruc.arrow_prop:
                if line_minus_arrow.find(arrow_) >= 0:
                    msg = self.error(line_no, "Extra arrow")
                    raise SyntaxError(msg)

            educts = _splitAtPlus(educts_str)
            products = _splitAtPlus(products_str)

            # Fix metabolite names (replace spaces with underscores)
            educts = self.fix_metab_names(educts, metabDict)
            products = self.fix_metab_names(products, metabDict)

            # Write line to output file
            outf.write(prefix + " + ".join(educts) + arrowStr +
                       " + ".join(products) + '\n')
        msg = ("Fixed %u reaction name(s) and %u metabolite name(s)." %
               (reaNamesFixed, len([key for key in metabDict
                                    if key != metabDict[key]])))
        self.messages.append((Verbosity.INFO, msg))

    def fix_metab_names(self, metab_strings, metabDict):
        """ process a list of metabolites; called internally by fixfileByHandles

        This function replaces whitespace in metabolite names with underscores
        and converts metabolite names to all lowercase. It tries to interpret
        the first "word" of a given metabolite name containing spaces as the
        stoichiometric coefficient.

        Keyword arguments:

        metab_strings -- list of metabolites given as strings in the form
                         "<coefficient> <metabolite>" (where metabolite may
                         still contain whitespace characters)
        metabDict     -- dictionary {name : fixed name}

        Returns:

        names_fixed   -- list of fixed "<coefficient> <metabolite>" strings,
                         where metabolite does not contain any whitespace
                         characters
        """
        newNames = set(metabDict.values())
        names_fixed = []
        for metab_string in metab_strings:
            if not metab_string:
                continue

            entry = metab_string.lower().split(None, 1)
            if len(entry) > 1:
                # Try to interpret first part as coefficient
                try:
                    float(entry[0])
                    names_fixed.append(entry[0] + ' ' +
                                       self.getFixedMetabName(entry[1], metabDict, newNames))
                except ValueError:
                    names_fixed.append(self.getFixedMetabName(
                        metab_string.lower(), metabDict, newNames))
            else:
                names_fixed.append(entry[0])
        return names_fixed

    def getFixedMetabName(self, name, metabDict, newNames):
        """ return fixed metabolite name; called internally by fix_metab_names()

        Keyword arguments:

        name        -- metabolite name to fix
        metabDict   -- dictionary {name : fixed name}
        newNames    -- set of metabDict.values()
        """
        if name in metabDict:
            return metabDict[name]
        fixed_name = name.replace(' ', '_')
        while fixed_name in newNames:
            fixed_name += '_'
        metabDict[name] = fixed_name
        newNames.add(fixed_name)
        if name != fixed_name:
            msg = "Metabolite '%s' --> '%s'" % (name, fixed_name)
            self.messages.append((Verbosity.DEBUG, msg))
        return fixed_name
