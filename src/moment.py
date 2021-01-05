#!/usr/bin/env/python

""" MOMENT: Modular Metabolic Engineering Tool

This script provides a modular metabolic engineering approach
using MOMA analysis. It iterates through a set of given reactions
and choose the best one for product yield enhancement based on
MOMA flux calculation. Afterwards the reaction is added to the
set of optimal reactions and the calculation goes on.

Runtime for a set of seven reactions: ~7 hours

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
from __future__ import print_function
from __future__ import division

from builtins import input
from builtins import zip
from builtins import str
from builtins import object
from past.utils import old_div
from metano.metabolicmodel import MetabolicModel
from metano.paramparser import ParamParser
from metano.defines import FbaParam
from metano.fba import FbAnalyzer
from metano.moma import MomaAnalyzer
from metano.metabolicflux import MetabolicFlux
from metano.moma import main
from multiprocessing import Pool
import os, time, glob, optparse, sys, json, errno


version_number = "alpha_1.3"

# ============================ SETTING PARAMETERS =============================== #

def unwrap_self_metabolic_engineering(arg, **kwarg):
    return runMoment.metabolic_engineering(*arg, **kwarg)


class OptionParser(optparse.OptionParser):
    def check_required(self, opt):
        option = self.get_option(opt)
        # Assumes the options default is set to None!
        if getattr(self.values, option.dest) is None:
            self.error("%s option not supplied" % option)
            sys.exit()




class runMoment(object):

    beginTime = time.time()

    def __init__(self, reactionFile, scenarioFile, outputFile, config, debugging=False):
        """
        initialize all required values for MOMENT
        """
        self.config = config
       # define the path to the reaction file:
        self.reactionFile = reactionFile
        # define the path to the scenario file:
        self.scenarioFile = scenarioFile
        # define the path to the output file:
        self.outputFile = outputFile
        # define the path to the output file for all best reactions:
        self.maxkeyOutputFile = outputFile[0:-4] + "_maxkey.txt"
        # define the path to the result folder:
        self.resultFolder = config["resultFolder"]
        # toggle debug mode, if true FBA is performed instead of MOMA:
        self.debugging = debugging
        self.numberMaxThreads = config["numberMaxThreads"]
        # define how many iteration cycles should be calculated:
        self.numberMetabolicInterventions = config["numberMetabolicInterventions"]
        # define how many reactions are considered as best for next cycle:
        #self.numberNextTurnReactions = config["numberNextTurnReactions"]
        # define the flux to be optimized and the wild type target flux:
        self.optimizationFlux = config["targetFlux"]
        #
        self.transporterPrefix = config["transporterPrefix"]
        #
        self.ignoreReactionPrefix = config["ignoreReactionPrefix"]
        self.treshold = config["minDiffTreshold"]
        # Set of reactions that are set as already optimized:
        self.setOfInterventions = config["setOfInterventions"]

        self.wtFbaFlux = self.Fba_solution(self.scenarioFile)


        # clear all previous parameters and result files and create if not exist:
        try:
            os.makedirs(self.resultFolder)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            for element in glob.glob(self.resultFolder + "*.txt"):
                os.remove(element)
        try:
            os.makedirs("insce/")
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
            for element in glob.glob("insce/*.txt"):
                os.remove(element)

        # read out the scenario file:
        cfg = open(self.scenarioFile, "r")
        self.insce = ""
        for line in cfg:
            if line[0] not in ("\n"):
                self.insce += line
        cfg.close()


    def Fba_solution(self, sceFile):
        """
        Precalculate the wild type solution for MOMA analysis
        """
        # read out the reactionsfile and parse the metabolic model:
        model = MetabolicModel()
        model.addReactionsFromFile(self.reactionFile)
        pparser = ParamParser()

        # set constraints based on modified scenario file:
        maxmin, objStr, solver, numIter, lb , ub = pparser.parse(sceFile)
        fbaParams = FbaParam(solver, maxmin, objStr, numIter)
        fbaParams.setLinConstraints(pparser.lin_constraints)

        scenarioModel = MetabolicModel()
        scenarioModel.addReactions(model.reactions)
        scenarioModel.setFiniteBounds(lb, ub, False)

        # perform FBA on the model with given constraints:
        fba = FbAnalyzer()
        obj_value, flux, _ = fba.runOnModel(scenarioModel, fbaParams, rmDeadEnds=True)
        return flux


    def identify_candidates(self):
        """
        identify all reaction candidates automatically^or load from file
        """

        # load input file which includes the set of reactions to be modified with
        # their values, seperated by colons:
        self.reactionsDct = {}
        self.setOfCandidates = []
        if len(self.config["setOfReactionCandidates"]) > 0:
            for item in self.config["setOfReactionCandidates"]:
                self.setOfCandidates.append(str(item.split(":")[0].strip()))
                self.reactionsDct[str(item.split(":")[0].strip())] = float(item.split(":")[-1].strip())
        else:
        # 1. perform FBA with optimizationFlux:
            subst =  "OBJ MAX " + self.optimizationFlux + " \n"

            # edit scenario file and write changes to new file:
            with open("sce.txt") as cfg:
                for line in cfg:
                    if line[0] == "#":
                        continue
                    if "OBJ " not in line:
                        subst += line
            scenario = open("insce/optimizationTarget.txt", "w")
            scenario.write(subst)
            scenario.close()
            self.optFbaFlux = self.Fba_solution("insce/optimizationTarget.txt")

            # 2. compare targetFBA with WT FBA:
            comparisonDict = {}
            for rea in self.wtFbaFlux:
                a = self.optFbaFlux[rea]
                b = self.wtFbaFlux[rea]
                diff = max(a, b) - min(a, b)
                comparisonDict[rea] = diff

            # 3. remove all specific reactions (diff < Th; transporter; production pathway)
            delete = []
            for rea in comparisonDict:
                if comparisonDict[rea] < self.treshold:
                    delete.append(rea)
                elif self.transporterPrefix in rea:
                    delete.append(rea)
                else:
                    for prefix in self.ignoreReactionPrefix:
                        if prefix in rea:
                            delete.append(rea)

            for rea in delete:
                del comparisonDict[rea]
            # check how many candidate reactions were calculated and
            # give the possibility to abort and adjust treshold
            if len(comparisonDict) > 30 or len(comparisonDict) < 5:
                print("""WARNING!! %d candidate reactions were calculated !\n
                Try to adjust the treshold in the config file!""" % len(comparisonDict))
                answer = input("Do you want to continue? (y/n)\n")
                if answer == "n":
                    exit()

            # 4. export values into dict:
            for rea in comparisonDict:
                self.setOfCandidates.append(rea)
                self.reactionsDct[rea] = float(self.optFbaFlux[rea])

        for element in self.setOfInterventions:
            self.setOfCandidates.remove(element)
        print("Candidate Reactions:", self.setOfCandidates)



    def metabolic_engineering(self,argument):
        """
        calculate the MOMA result for target flux while setting constraints
        to all optimization reactions.
        """
        new_rea, setOfInterventions_list = argument
        name = str(len(setOfInterventions_list)) + "_" + new_rea

        if new_rea == []:
            pass
        else:
            # modify the scenario file by setting upper and lower bounds
            # for each manipulated reaction:
            subst =  "LB " + new_rea + " " + str(self.reactionsDct[new_rea]) + " \n"
            subst += "UB " + new_rea + " " + str(self.reactionsDct[new_rea]) + " \n\n"
            for element in setOfInterventions_list:
                subst += "LB " + element + " " + str(self.reactionsDct[element]) + " \n"
                subst += "UB " + element + " " + str(self.reactionsDct[element]) + " \n"

            # write changes to scenario file:
            cfg = open("insce/insce_" + name + ".txt", "w")
            cfg.write(self.insce + "\n\n")
            cfg.write(subst)
            cfg.close()

            # read out the reactionsfile and parse the metabolic model:
            model = MetabolicModel()
            model.addReactionsFromFile(self.reactionFile)
            pparser = ParamParser()

            # set constraints based on modified scenario file:
            maxmin, objStr, solver, numIter, lb , ub = pparser.parse("insce/insce_" + name + ".txt")
            fbaParams = FbaParam(solver, maxmin, objStr, numIter)
            fbaParams.setLinConstraints(pparser.lin_constraints)
            scenarioModel = MetabolicModel()
            scenarioModel.addReactions(model.reactions)
            scenarioModel.setFiniteBounds(lb, ub, False)

            #save flux distribution of WT FBA into variable for usage during MOMA analysis:


            # during debugging, only FBA is performed:
            if self.debugging == True:
                solutionFlux = self.Fba_solution("sce.txt")
            else:
                # perform MOMA analysis and save the resulting flux
                # distribution in "solutionFlux":
                moma = MomaAnalyzer("default")
                dist, solutionFlux, status = moma.runOnModel(scenarioModel,
                    self.wtFbaFlux, fbaParams.linConstraints, numIter)[:3]


            # save the flux through optimization flux
            try:
                targetFlux = solutionFlux.fluxDict[self.optimizationFlux]
            except:
                targetFlux = 0
                print("Warning:", self.optimizationFlux, "could not be found!")

            # write result to file. This step is necessary when
            # using multiprocessing. In this way the cores dont
            # get in each others way:
            maxFluxFile = open(self.resultFolder + name + ".txt", "w")
            output = new_rea + " " + str(targetFlux) + "\n"
            maxFluxFile.write(output)

            # delete the scenario file:
            #os.remove("self.insce/self.insce_" + name + ".txt")
            return targetFlux


    def moment(self):

        # create empty output file:
        maxkey_output = ""
        text_file = open(self.outputFile, "w")
        text_file.write("Number;Reaction;Target_flux\n")
        text_file.close()

        # set number of cores to be used for multiprocessing:
        pool=Pool(processes=self.numberMaxThreads)

        # perform calculation of the iteration cycle:
        while len(self.setOfInterventions) < self.numberMetabolicInterventions:
            # set the calculation parameter and perform MOMA analysis:
            maxFluxDict = {}
            modify = []
            for entry in self.setOfCandidates:
                modify.append((entry, self.setOfInterventions))

            if self.debugging == True:
                pool=Pool(processes=1)
            pool.map(unwrap_self_metabolic_engineering, list(zip([self]*len(modify), modify)))
            #pool.map(self.metabolic_engineering, modify)
            # read out all results from primary files:
            for element in self.setOfCandidates:
                with open(self.resultFolder + str(len(self.setOfInterventions)) + "_" + element + ".txt") as f:
                    for line in f:
                        print(line)
                        maxFluxDict[line.split()[0]] = float(line.strip("\n").split()[-1])

            # write all results to output file:
            text_file = open(self.outputFile, "a")
            output = ""
            for key in maxFluxDict:
                output += str(len(self.setOfInterventions)) + ";" + key + ";" + str(maxFluxDict[key]) + "\n"
            text_file.write(output)
            text_file.close()
            print(output)

            # save the best result into self.setOfInterventions list and delete from new_rea:
            print(max(maxFluxDict, key=maxFluxDict.get))
            maxkey = max(maxFluxDict, key=maxFluxDict.get)
            print("***   ", maxkey, "=", maxFluxDict[maxkey], "   ***\n\n\n")
            self.setOfInterventions.append(maxkey)
            self.setOfCandidates.remove(maxkey)
            maxkey_output += maxkey + "; " + str(maxFluxDict[maxkey]) + "\n"

            print("Runtime for", len(self.setOfInterventions), "reaction(s):\n", \
                  round(old_div((time.time() - self.beginTime),60), 2), "minutes\n")


        # finally write the set of best reactions to output file:
        text_file = open(self.maxkeyOutputFile, "w")
        text_file.write(maxkey_output)
        text_file.close()
        print("\n\nMAXKEY OUTPUT:\n", maxkey_output)

        print("RUNTIME:\n", \
              round(time.time() - self.beginTime, 2), "seconds\n", \
              round(old_div((time.time() - self.beginTime),60), 2), "minutes\n", \
              round(old_div((time.time() - self.beginTime),60/60), 2), "hours\n")


def main():
    usage = "Usage: %prog [options]"
    version = "Modular Metabolic Engineering Test \nVersion %s" % version_number
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform Flux Balance Analysis on the network given "
                           "by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="scenarioFile",
                      help="use the given WT scenario FILE for Flux Balance "
                           "Analysis", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write output of MOMENT to FILE",
                      metavar="FILE")
    parser.add_option("-c", "--config", dest="configFile",
                      help="read all configurations from JSON FILE",
                      metavar="FILE")
    #parser.add_option("-t", "--target", action="store", type="string",
    #                  dest="targetFlux", help="target flux that should be maximized")
    parser.add_option("-D", "--debug", action="store_true",
                      dest="debugging", help="uses FBA instead of MOMA for shorter Runtime")
    #parser.add_option("-N", "--number", action="store", type="int",
    #                  dest="numberMetabolicInterventions", help="number of modified reactions")
    parser.set_defaults(debugging=False)

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-o")
    parser.check_required("-c")

    with open(options.configFile) as cf:
        cfDict = json.load(cf)

    myMoment = runMoment(options.reactionFile,
        options.scenarioFile,
        options.outputFile,
        cfDict,
        options.debugging)

    myMoment.identify_candidates()
    myMoment.moment()

if __name__ == "__main__":
    main()

