
'''
This module evaluates the objective function by setting each
metabolite as export objective function.
As result this file prints all metabolites that could not be
produced by the network.

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
'''
from __future__ import print_function

from builtins import object
import optparse
import os

from metano.metabolicmodel import MetabolicModel
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.fba import FbAnalyzer
from metano.defines import FbaParam, COPYRIGHT_VERSION_STRING


class OptionParser(optparse.OptionParser):
    """ Extension of OptionParser with function for checking required arguments
    """

    def check_required(self, opt):
        """ check if required options are given """
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            self.error("%s option not supplied" % option)


class EvalObjective(object):
    """ Evaluates objective function by testing each metabolite """
    def __init__(self, model, scenario):
        self.model = model
        self.scenario = scenario
        self.run_count = 0
        self.result = []

        objfunindex = model.reactionDict[model.biomass_name]
        metabolites = model.reactions[objfunindex].metabolites

        self.objfun_dict = {}
        for met in metabolites:
            self.objfun_dict[met.name] = met.coef

    def run(self):
        """ start evaluation process and save results """
        print ("Start evaluation of %s metabolites from "
               "the objective function." % len(self.objfun_dict))
        for metabolite in self.objfun_dict:
            pparser = ParamParser()
            maxmin, objstr, solver, numiter, lb , ub = pparser.parse(self.scenario)

            scenariomodel = MetabolicModel()
            scenariomodel.addReactions(self.model.reactions)
            scenariomodel.addReaction(("temp_obj_fun", ((1, metabolite),), "-->", ()))
            #objfunindex = scenariomodel.reactionDict["temp_obj_fun"]

            fbaparams = FbaParam(solver, True, "temp_obj_fun", numiter)
            fbaparams.setLinConstraints(pparser.lin_constraints)

            scenariomodel.setFiniteBounds(lb, ub, True)
            fba = FbAnalyzer(verbose=False)
            obj_value, solution, _ = fba.runOnModel(scenariomodel, fbaparams, rmDeadEnds=True)
            if not solution:
                self.result.append(metabolite)


    def print_results(self):
        """ print results from run """
        if not self.result:
            print("All metabolites from the objective functions are produced by the network!")
            return
        print("-----------------------------------------------------------------------")
        print("RESULT: the following metabolites could not be produced in the network:")
        for i in sorted(self.result):
            print("\t" + i)
        print("-----------------------------------------------------------------------")


def main():
    """ main function, parse command line and pass options to class """
    # 1. Parse command line

    usage = "Usage: %prog [options]"
    version = "multiple Flux balance analysis\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform multiple Flux Balance Analysis on the network given "
                           "by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="basic scenario file", metavar="FILE")

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    #parser.check_required("-f")

    # 1. parse reaction file
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

    # 2. evaluate objective function
    evalobj = EvalObjective(model, options.paramFile)
    evalobj.run()
    evalobj.print_results()



if __name__ == "__main__":
    main()
