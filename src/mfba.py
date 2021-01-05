
''' Performs batch computation of flux balance analysis

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
'''
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import optparse
import os
from multiprocessing import Pool

from metano.metabolicmodel import MetabolicModel
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.fba import FbAnalyzer
from metano.defines import FbaParam, COPYRIGHT_VERSION_STRING


def unwrap_self_run(arg, **kwarg):
    """ snipped for running multiprocessing within a class
    from Thomas Rueckstiess:
    http://www.rueckstiess.net/research/snippets/show/ca1d7d90
    """
    return MultipleFBA.perform_fba(*arg, **kwarg)


class OptionParser(optparse.OptionParser):
    """ Extension of OptionParser with function for checking required arguments
    """

    def check_required(self, opt):
        """ check required arguments """
        if type(opt) == list:
            for i in opt:
                option = self.get_option(i)
                if getattr(self.values, option.dest) is not None:
                    print(i, getattr(self.values, option.dest))
                    return
            self.error(
                "you must provide at least on of these options: %s" % ", ".join(opt))

        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            self.error("%s option not supplied" % option)


class MultipleFBA(object):
    """ Class for batch computation of FBA """

    def __init__(self, model, scenario, constraints, fba_file_prefix,
                 reactions_of_interest, verb=False):
        """ init """
        self.model = model
        self.verb = verb
        self.insce = []
        self.rates = {}
        self.run_count = 0
        self.fba_file_prefix = fba_file_prefix
        self.result = []
        self.roi = reactions_of_interest

        if constraints[0].split("\t")[0] != "NAME":
            print ("Error in file %s:"
                   "First line must begin with NAME flag!"
                   "See Metano documentation for further information")
            exit()
        self.names = [i.strip() for i in constraints[0].split("\t")[1:]]
        duplicates = set([x for x in self.names if self.names.count(x) > 1])
        if duplicates:
            print("Error: Found duplicate scenario name(s) in constraint matrix:")
            print("---", ", ".join(duplicates), "---")
            print("Scenario names MUST be unique!")
            exit()
        del constraints[0]

        for line in constraints:
            line = [i.strip() for i in line.split("\t")]
            if line[0] not in ("\t", "\n", "", " "):
                self.rates[line[0]] = [i if i else "0.0" for i in line[1:]]

        try:
            self.run_count = len(line)-1
        except UnboundLocalError:
            print("Error: Constraint matrix is empty!")
            exit()

        for line in scenario:
            if line[0] not in "\n":
                self.insce.append(line)

    def perform_fba(self, run):
        """ performs multiple fba """
        insce_temp = self.insce[:]
        for constraint in self.rates:
            insce_temp.append(
                "\t".join([constraint, self.rates[constraint][run]]))

        pparser = ParamParser()

        maxmin, obj_str, solver, num_iter, lb, ub = pparser.parseByHandle(
            insce_temp)

        fba_params = FbaParam(solver, maxmin, obj_str, num_iter)
        fba_params.setLinConstraints(pparser.lin_constraints)
        scenario_model = MetabolicModel()
        scenario_model.addReactions(self.model.reactions)
        scenario_model.setFiniteBounds(lb, ub, True)
        fba = FbAnalyzer()
        obj_value, solution, _ = fba.runOnModel(
            scenario_model, fba_params, rmDeadEnds=True)

        result = [obj_value]
        messages = []
        for rea in self.roi:
            try:
                result.append(solution.fluxDict[rea])
            except KeyError:
                messages.append(
                    "WARNING: %r could not be found in reactions" % rea)
                result.append(0.0)
                # self.roi.remove(rea)
        if messages and len(solution) != 0:
            print("-" * max([len(i) for i in messages]))
            print("\n".join(messages))
            print("-" * max([len(i) for i in messages]))

        if self.verb:
            print(self.names[run] + ":")
        if len(solution) != 0:
            if self.verb:
                print("Value of obj function:", obj_value)
            if self.fba_file_prefix:
                try:
                    solution.writeToFile(
                        self.fba_file_prefix + str(self.names[run]) + ".txt")
                except IOError as strerror:
                    print ("Unable to write to file %s:" %
                           os.path.basename(self.fba_file_prefix + str(self.names[run]) + ".txt"))
                    print(strerror)
                    exit()
        elif self.verb:
            print("No output written.")
        return (self.names[run],) + tuple(result)

    def run(self):
        """ performs one run at a time """
        result = []
        for run in range(self.run_count):
            res = self.perform_fba(run)
            result.append(res)
        return result

    def run_pool(self, cores):
        """ use multiprocessing for batch computation """
        pool = Pool(processes=cores)
        result = pool.map(unwrap_self_run, list(
            zip([self]*self.run_count, list(range(self.run_count)))))
        return result

    def write_solution_to_file(self, filename, solution):
        """ write solution to file """
        write_solution = [
            "\t".join(["Scenario_name", "objective"] + [i for i in self.roi])]
        for line in solution:
            write_solution.append("\t".join([str(i) for i in line]))
        with open(filename, "w") as f:
            f.write("\n".join(write_solution))


def main():
    """ parse the command line, create the mFBA class and call anylsis functions """

    usage = "Usage: %prog [options]"
    version = "multiple Flux balance analysis\n" + COPYRIGHT_VERSION_STRING
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--reactions", dest="reactionFile",
                      help="perform multiple Flux Balance Analysis on the network given "
                           "by the reaction FILE", metavar="FILE")
    parser.add_option("-p", "--parameters", dest="paramFile",
                      help="basic scenario file", metavar="FILE")
    parser.add_option("-m", "--matrix", dest="matrixFile",
                      help="matrix that contains all constraints", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="write CSV result to file, tab-seperated", metavar="FILE")
    parser.add_option("-f", "--solution-files", dest="fluxFilePrefix",
                      help="write FBA solutions to files marked "
                      "with FILE-PREFIX (may contain path)",
                      metavar="FILE-PREFIX")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verb", help="print the obj function of each FBA solution")
    parser.add_option("-c", "--cores", type="int",
                      dest="cores", help="defines how many CPU cores to use; "
                      "by default multiprocessing is disabled")
    parser.add_option("-l", "--list-of-reactions", dest="reactionsOfInterest",
                      help="a LIST of additional reactions of interest, seperated by comma; "
                      "if nothing is provided, only the obj function is written "
                      "into the result file")
    parser.set_defaults(verb=False, cores=0, reactionsOfInterest="")

    options, _ = parser.parse_args()
    parser.check_required("-r")
    parser.check_required("-p")
    parser.check_required("-m")
    parser.check_required(["-f", "-o"])

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

    if options.reactionsOfInterest:
        reactions_of_interest = [i.strip()
                                 for i in options.reactionsOfInterest.split(",")]
    else:
        reactions_of_interest = []

    matrix = open(options.matrixFile, "r")
    matrix = matrix.readlines()
    scenario = open(options.paramFile, "r")

    mfba = MultipleFBA(model, scenario, matrix, options.fluxFilePrefix,
                       reactions_of_interest, options.verb)

    if options.cores == 0:
        result = mfba.run()

    else:
        result = mfba.run_pool(options.cores)

    if options.outputFile:
        print("write ouput to file")
        mfba.write_solution_to_file(options.outputFile, result)


if __name__ == "__main__":
    main()
