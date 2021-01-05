#!/usr/bin/env python

"""
This script predicts the directions of reactions based on the upper and lower
bounds for the Gibbs free energy changes computed by getReaEnthalpies(). It
checks whether the given reaction and scenario file conform to the predicted
reaction directions and adds appropriate constraints to the given scenario file
to make reactions irreversible.


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

from metano.defines import COPYRIGHT_VERSION_STRING
from metano.fba import OptionParser
from metano.optparse import OptionGroup
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.metabolicmodel import MetabolicModel
from metano.compute_rea_enthalpies import getReaEnthalpies, readReaEnthalpiesFromFile
import os


def main():

    # Parse command line
    usage = ("Usage: %prog <delta-G-file> [options]\n\n"
             "Gibbs free energy changes are taken from the <delta-g-file> if "
             "given.\nElse, they are computed from the data specified via the "
             "-r/-c/-d/-s/-t set of\noptions.")
    version = "%prog\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    standardOptGroup = OptionGroup(parser, "Standard parameters")
    standardOptGroup.add_option("-p", "--parameters", dest="paramFile",
                                help="check the given scenario FILE",
                                metavar="FILE")
    standardOptGroup.add_option("-o", "--output", dest="outputFile",
                                help="write modified scenario file to FILE "
                                "(must not be the same as scenario file)",
                                metavar="FILE")
    standardOptGroup.add_option("-e", "--epsilon", dest="epsilon",
                                help="set THRESHOLD for recognizing Gibbs free "
                                "energy change as clearly positive (default 0)",
                                metavar="THRESHOLD")
    parser.add_option_group(standardOptGroup)

    computeOptGroup = OptionGroup(parser, "Parameters for computation of "
                                  "Gibbs free energies", "These are only needed"
                                  " if delta-G values are not read from file.")
    computeOptGroup.add_option("-r", "--reactions", dest="reactionFile",
                               help="use reactions from reaction FILE",
                               metavar="FILE")
    computeOptGroup.add_option("-c", "--concentrations",
                               dest="concentrationFile", help="use "
                               "concentrations from FILE", metavar="FILE")
    computeOptGroup.add_option("-d", "--thermodynamics", dest="thermodynFile",
                               help="use thermodynamic data from FILE",
                               metavar="FILE")
    computeOptGroup.add_option("-s", "--synonyms", dest="synonymFile",
                               help="use metabolite synonyms from FILE",
                               metavar="FILE")
    computeOptGroup.add_option("-t", "--temperature", dest="temperature",
                               help="set temperature to VALUE (in degrees "
                               "Celsius)", metavar="VALUE")
    parser.add_option_group(computeOptGroup)
    parser.set_defaults(epsilon='0.')

    options, args = parser.parse_args()
    parser.check_required("-p")
    parser.check_required("-o")
    parser.check_required("-t")

    try:
        epsilon = float(options.epsilon)
    except ValueError:
        print ("Error: Invalid floating point value for epsilon (%s)" %
               options.epsilon)
        exit()
    if (os.path.exists(options.outputFile) and
            os.path.samefile(options.outputFile, options.paramFile)):
        print ("Error: Input and output scenario files are the same (%s)" %
               options.paramFile)
        exit()

    if epsilon < 0.:
        print("Warning: epsilon < 0. Using default value of 0 instead.")
        epsilon = 0.

    if len(args) > 0:
        gibbsR = readReaEnthalpiesFromFile(args[0])
    else:
        print ("\nInfo: No file with Gibbs free energies given. Launching "
               "computation of values.\n")
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
        gibbsR = getReaEnthalpies(options.concentrationFile,
                                  options.thermodynFile, options.synonymFile,
                                  options.reactionFile, temperature)

    # Parse scenario file
    pparser = ParamParser()
    try:
        # Parse file
        maxmin, obj_name, solver, numiter, lb, ub =\
            pparser.parse(options.paramFile)
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

    if options.reactionFile:
        # Parse reaction file
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

        # Set flux bounds in model
        model_messages = []
        try:
            model.setFiniteBounds(lb, ub, True, model_messages)
        except ValueError as strerror:
            print(strerror)
            exit()
        reactions = model.getReactionNames()
        # Include implicit flux bounds from reaction file
        p_lb, p_ub = lb, ub
        lb, ub = model.getFiniteBounds()

        # Show warning and info messages of parsers
        msgs = (rparser.getMessages() + pparser.getMessages() +
                [x[1] for x in model_messages])
        if msgs:
            print('\n'+'\n'.join(msgs))
    else:
        print ("Warning: No reaction file given. Implicit bounds defined by "
               "reaction\nreversibilities cannot be taken into account.")
        p_lb, p_ub = lb, ub
        reactions = sorted(
            set(list(gibbsR.keys())+list(lb.keys())+list(ub.keys())))

    # Copy optimization parameters
    param_head = ("OBJ %s %s\nSOLVER %s\n" % ({True: "MAX", False: "MIN"}[maxmin],
                                              obj_name, solver))
    if numiter >= 0:
        param_head += "NUM_ITER %u\n" % numiter
    param_head += '\n'

    try:
        with open(options.outputFile, 'w') as f:
            f.write(param_head)
            for rea in reactions:
                if rea in gibbsR:
                    mini, maxi = gibbsR[rea]
                    if maxi < -epsilon:
                        # Reaction is irreversible and proceeds from
                        # left to right (LB 0)
                        if rea in ub and ub[rea] <= 0.:
                            print ("Warning: Gibbs free energy predicts "
                                   "LB %s 0,\nbut reaction is defined as UB %g "
                                   "(skipped)" % (rea, ub[rea]))
                            if rea in p_lb:
                                f.write("LB %s %g\n" % (rea, p_lb[rea]))
                            if rea in p_ub:
                                f.write("UB %s %g # Conflict (delta_G_max = %g)"
                                        "\n" % (rea, p_ub[rea], maxi))
                        elif rea in lb and lb[rea] >= 0.:
                            if rea in p_lb:
                                f.write("LB %s %g # Direction confirmed by "
                                        "delta_G_max = %g\n" % (rea, p_lb[rea],
                                                                maxi))
                            if rea in p_ub:
                                f.write("UB %s %g\n" % (rea, p_ub[rea]))
                        else:
                            # (rea not in ub or ub[rea] > 0) and (rea not in lb
                            # or lb[rea] < 0)
                            # Add LB 0 line
                            f.write("LB %s %g # Assigned because delta_G_max = "
                                    "%g\n" % (rea, 0., maxi))
                            if rea in p_ub:
                                f.write("UB %s %g\n" % (rea, p_ub[rea]))
                    elif mini > epsilon:
                        # Reaction is irreversible and proceeds from
                        # right to left (UB 0)
                        if rea in lb and lb[rea] >= 0.:
                            print ("Warning: Gibbs free energy predicts "
                                   "UB %s 0,\nbut reaction is defined as LB %g "
                                   "(skipped)" % (rea, lb[rea]))
                            if rea in p_lb:
                                f.write("LB %s %g # Conflict (delta_G_min = %g)"
                                        "\n" % (rea, p_lb[rea], mini))
                            if rea in p_ub:
                                f.write("UB %s %g\n" % (rea, p_ub[rea]))
                        elif rea in ub and ub[rea] <= 0.:
                            if rea in p_lb:
                                f.write("LB %s %g\n" % (rea, p_lb[rea]))
                            if rea in p_ub:
                                f.write("UB %s %g # Direction confirmed by "
                                        "delta_G_min = %g\n" % (rea, p_lb[rea],
                                                                mini))
                        else:
                            # (rea not in ub or ub[rea] > 0) and (rea not in lb
                            # or lb[rea] < 0)
                            if rea in p_lb:
                                f.write("LB %s %g\n" % (rea, p_lb[rea]))
                            # Add UB 0 line
                            f.write("UB %s %g # Assigned because delta_G_min = "
                                    "%g\n" % (rea, 0., mini))
                    else:
                        # Just copy parameters
                        if rea in p_lb:
                            f.write("LB %s %g\n" % (rea, p_lb[rea]))
                        if rea in p_ub:
                            f.write("UB %s %g\n" % (rea, p_ub[rea]))
                else:
                    # Just copy parameters
                    if rea in p_lb:
                        f.write("LB %s %g\n" % (rea, p_lb[rea]))
                    if rea in p_ub:
                        f.write("UB %s %g\n" % (rea, p_ub[rea]))
    except IOError as strerror:
        print ("Unable to write to file %s:" %
               os.path.basename(options.outputFile))
        print(strerror)
        exit()


if __name__ == "__main__":
    main()
