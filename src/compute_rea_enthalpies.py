#!/usr/bin/env python

"""
This script computes Gibbs free reaction energies from the Gibbs energies of
formation and the concentrations of the involved metabolites. It uses two values
for each concentration, an upper and a lower bound for the physiological levels
of each metabolite, and returns an upper and a lower bound for the Gibbs free
energy change for each reaction based on the concentration bounds.

It takes a file with the concentrations, a file with the Gibbs free energies of
formation, a synonym table (also from a file), a metano reaction file, and the
temperature as parameters.
The output is a file which lists the computed minimum and maximum values for the
Gibbs free energy change of each reaction.


This file is part of metano.
Copyright (C) 2010-2019 Alexander Riemer, Julia Helmecke
Braunschweig University of Technology,
Dept. of Bioinformatics and Biochemistry

Michael Stelzer and Rene Rex contributed to this file.

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

import csv
import re
import os
from sys import stdout
from collections import defaultdict
from math import log
from metano.defines import COPYRIGHT_VERSION_STRING
from metano.fba import OptionParser
from metano.reactionparser import ReactionParser
from metano.metabolicmodel import MetabolicModel


R = 8.314472           # J mol-1 K-1, universal gas constant
KelvinOffset = 273.15  # offset for conversion from Celsius to Kelvin


def _replaceChars(s):
    return re.sub(" |:", "_", s)


def getReaEnthalpies(concentrationFile, thermodynFile, synonymFile,
                     reactionFile, temperatureCelsius):

    # result dict {reaction name : (delta G min, delta G max)}
    result = {}

    # Read synonyms from file into dictionary {synonym : preferred name}
    with open(synonymFile) as f:
        synonymReader = csv.DictReader(f, delimiter='\t')

        synonyms = {}
        for row in synonymReader:
            preferredName = _replaceChars(row[synonymReader.fieldnames[0]])
            synonyms[preferredName] = preferredName
            synonyms[_replaceChars(row[synonymReader.fieldnames[1]])] =\
                preferredName
            brendaSyn = row[synonymReader.fieldnames[3]]
            keggSyn = row[synonymReader.fieldnames[4]]
            metaSyn = row[synonymReader.fieldnames[5]]
            if brendaSyn != '-' and brendaSyn != '':
                for syn in brendaSyn.split(", "):
                    synonyms[_replaceChars(syn)] = preferredName
            if keggSyn != '-' and keggSyn != '':
                for syn in keggSyn.split("##"):
                    synonyms[_replaceChars(syn.split("::")[1])] = preferredName
            if metaSyn != '-' and metaSyn != '':
                for syn in metaSyn.split("##"):
                    synonyms[_replaceChars(syn.split("::")[1])] = preferredName

    # Read concentrations from file into dict {linear expression : (min, max)}
    with open(concentrationFile) as f:
        concentrationReader = csv.DictReader(f, delimiter=';')

        concentrations = {}
        for row in concentrationReader:

            name = _replaceChars(row[concentrationReader.fieldnames[1]])
            name = synonyms.get(name, name)

            concentrations[name] = (float(row[concentrationReader.fieldnames[2]]),
                                    float(row[concentrationReader.fieldnames[3]]))

    # Read thermodynamic data from file into a
    # dict {metabolite : list of 4-tuples}
    with open(thermodynFile) as f:
        thermodynReader = csv.reader(f, delimiter=';')
        thermodyn = defaultdict(list)
        currentName = ""
        for row in thermodynReader:
            if row[0]:
                currentName = re.sub(r" \((aq|ion|red|ox|tot)\)", "", row[0])
                currentName = _replaceChars(currentName)
                currentName = synonyms.get(currentName, currentName)
            else:
                thermodyn[currentName].append(tuple(row[1:5]))
        # We don't want the defaultdict behavior anymore after this point
        thermodyn = dict(thermodyn)

    # Parse reaction file
    rparser = ReactionParser()
    model = MetabolicModel()
    try:
        model.addReactionsFromFile(reactionFile, rparser)
    except IOError as strerror:
        print ("An error occurred while trying to read file %s:" %
               os.path.basename(reactionFile))
        print(strerror)
        exit()
    except SyntaxError as strerror:
        print ("Error in reaction file %s:" %
               os.path.basename(reactionFile))
        print(strerror)
        exit()

    temperatureKelvin = temperatureCelsius + KelvinOffset

    for reaction in model:

        DeltaG0 = 0.
        logTermMax = 0.
        logTermMin = 0.

        success = True

        for coef, name in reaction:
            name = synonyms.get(name, name)
            try:
                DeltaG0 += coef * float(thermodyn[name][0][0])
                cMin, cMax = concentrations[name]

                if coef > 0:
                    logTermMax += coef * log(cMax)
                    logTermMin += coef * log(cMin)
                else:
                    logTermMax += coef * log(cMin)
                    logTermMin += coef * log(cMax)

            except KeyError:
                success = False
                break

        if success:
            DeltaGmin = DeltaG0 + R*temperatureKelvin*logTermMin
            DeltaGmax = DeltaG0 + R*temperatureKelvin*logTermMax
            result[reaction.name] = DeltaGmin, DeltaGmax

    return result


def readReaEnthalpiesFromFile(filename):
    with open(filename) as f:
        return readReaEnthalpiesFromFileHandle(f)


def readReaEnthalpiesFromFileHandle(f):
    result = {}
    line_no = 0
    for line in f:
        line_no += 1
        if line_no == 1:
            # Skip first line (table header)
            continue
        key, mini, maxi = line.split()
        result[key] = float(mini), float(maxi)
    return result


def writeOutputToFileHandle(f, solution):
    if len(solution) == 0:
        return  # Nothing to do
    # Get length of longest reaction name (for formatting)
    maxlen = len(max(solution, key=len))
    maxlen = max(maxlen, len("REACTION"))
    # Write table head
    f.write("REACTION".ljust(maxlen)+" "+"MIN_DELTA_G".rjust(12)+" " +
            "MAX_DELTA_G".rjust(12)+"\n")
    # Write Gibbs free energy bounds for each reaction
    for rea in sorted(solution):
        mini, maxi = solution[rea]
        f.write(rea.ljust(maxlen) + " " +
                ("% .6g" % mini).rjust(12) + " " +
                ("% .6g" % maxi).rjust(12) + "\n")


def main():
    # Parse command line
    usage = "Usage: %prog [options]"
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="use reactions from reaction FILE",
                      metavar="FILE")
    parser.add_option("-c", "--concentrations", dest="concentrationFile",
                      help="use concentrations from FILE",
                      metavar="FILE")
    parser.add_option("-d", "--thermodynamics", dest="thermodynFile",
                      help="use thermodynamic data from FILE",
                      metavar="FILE")
    parser.add_option("-s", "--synonyms", dest="synonymFile",
                      help="use metabolite synonyms from FILE",
                      metavar="FILE")
    parser.add_option("-t", "--temperature", dest="temperature",
                      help="set temperature to VALUE (in degrees Celsius)",
                      metavar="VALUE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write output (Gibbs free energies) to FILE",
                      metavar="FILE")

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-c")
    parser.check_required("-d")
    parser.check_required("-s")
    parser.check_required("-t")

    try:
        temperature = float(options.temperature)
    except ValueError:
        print ("Error: Invalid floating point value for temperature (%s)" %
               options.temperature)
        exit()

    # Compute Gibbs free energies from the given data
    gibbsR = getReaEnthalpies(options.concentrationFile, options.thermodynFile,
                              options.synonymFile, options.reactionFile,
                              temperature)

    # Write output either to a given file or to the console
    if options.outputFile:
        try:
            with open(options.outputFile, 'w') as f:
                writeOutputToFileHandle(f, gibbsR)
        except IOError as strerror:
            print ("Unable to write to file %s:" %
                   os.path.basename(options.outputFile))
            print(strerror)
            writeOutputToFileHandle(stdout, gibbsR)
    else:
        writeOutputToFileHandle(stdout, gibbsR)


if __name__ == "__main__":
    main()
