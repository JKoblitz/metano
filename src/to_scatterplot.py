#!/usr/bin/env python

""" create scatterplot from scenario file

This script will create a scatterplot from given flux distributions or
scenario files.

It parses a reaction file and a scenario file, formulates the
corresponding linear, quadratic, or other nonlinear optimization problem,
and calls an appropriate solver.


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

from builtins import zip
from builtins import range
from builtins import object
from past.utils import old_div
import os
import random
import re
import argparse
import configparser
import matplotlib.pyplot as plt
from metano.metabolicflux import MetabolicFlux
from metano.fba import FbAnalyzer
from metano.fva import FvAnalyzer
from metano.mfm import MfmAnalyzer
from metano.metabolicmodel import MetabolicModel
from metano.reactionparser import ReactionParser
from metano.paramparser import ParamParser
from metano.defines import FbaParam, COPYRIGHT_VERSION_STRING

try:
    from adjustText import adjust_text
    ADJ_TEXT_AV = True
except ImportError:
    ADJ_TEXT_AV = False
    print("""Warning: The 'adjustText' package is not installed.
             It is highly recommended to use adjustText when annotations are shown,
             since the annotation placement is significantly improved.\n
             To install adjustText, please type:\n
             (sudo) pip install adjustText""")

ALLOWED_METHODS = ["fba", "fva", "mfm", "mva", "moma"]


def abs_diff(a_vec, b_vec):
    """ calculate the absolute difference between two vectors """
    if max(a_vec) > max(b_vec):
        return min(a_vec) - max(b_vec)
    return min(b_vec) - max(a_vec)


def abbreviate_all_in_list(str_list, max_len=20, ids=False, unique=False):
    """ abbreviate all strings in given list """
    result = []
    letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    special_chars = ('_', '-', '"', "'", ',')
    lower_vowels = ('a', 'e', 'i', 'o', 'u')
    upper_vowels = ('A', 'E', 'I', 'O', 'U')

    for item in str_list:
        found = False
        if ids:
            # use database ids if available
            regex_id = {"brenda": re.compile(r"(B[S|R]\d+)"),
                        "metacyc": re.compile(r"(RXN[A-Z\-\d]+|[A-Z\.\-\d\+]+-RXN)"),
                        "kegg": re.compile(r"(R\d+)")}

            for database in regex_id:
                if regex_id[database].findall(item):
                    found = True
                    result.append(regex_id[database].findall(item)[0])
                    break
        # if no id found or id search not applicable
        if not found:
            item = item.strip()
            if len(item) <= max_len:
                if item in result:
                    item = item[:-1] + random.choice(letters)
                else:
                    result.append(item)
                    continue
            while len(item) > max_len:
                # always start to remove items from the end of the string
                # first delete all lower vowels
                if any(l in item for l in lower_vowels):
                    index = max([item.rfind(l) for l in lower_vowels])
                    item = item[:index] + item[index+1:]
                # then delete all upper vowels
                elif any(l in item for l in upper_vowels):
                    index = max([item.rfind(l) for l in upper_vowels])
                    item = item[:index] + item[index+1:]
                # next delete all special characters
                elif any(l in item for l in special_chars):
                    index = max([item.rfind(l) for l in special_chars])
                    item = item[:index] + item[index+1:]
                else:
                    # at the end cut the remaining elements from max_len
                    item = item[0:max_len]
            # if each item should be unique in list:
            if unique:
                j = 0
                while item in result:
                    j += 1
                    item = item[:-1] + random.choice(letters)
                    if j > 20:
                        break
            result.append(item)
    return result


class Scatterplot(object):
    """ Class for the generation of a scatterplot"""

    def __init__(self, method, filenames, model=None):
        super(Scatterplot, self).__init__()

        self.method = method.lower()
        if self.method not in ALLOWED_METHODS:
            print ("Error: Method %r is unknown. Please choose one of the following "
                   "methods: FBA, FVA, MFM, MVA (only for sce)" % (self.method))
            exit()

        self.filenames = filenames
        if len(self.filenames) != 2:
            print("Error: -f option needs exactly 2 solution or "
                  "scenario files to create a scatterplot.")
            exit()

        self.model = model
        self.result = {}

        self.style = {
            "label_length": "15",
            "ids_in_names": "False",
            "unique_names": "False",
            "min_value": "0",
            "insignificant_color": "#545454",
            "significant_color": "#BE1E3C",  # TU BS red
            "xlabel": self.filenames[0],
            "ylabel": self.filenames[1],
            "title": "",
            "show_errorbar": "significant",  # significant, all, none
            "show_labels": "significant",  # significant, all, none
            "arrowcolor": "#BE1E3C",
            "labelcolor": "#BE1E3C",
            "significance_threshold": "0.1",
            "fontsize": "10",
            "figsize": "10 8",
            "plot_grid": "False",
            "arrow_line_width": "1",
        }

        if self.is_sce():
            if not self.model:
                print("Error: a reaction file is needed to perform "
                      "analysis from given scenario files!")
                exit()
            self.do_analysis()
        else:
            self.read_solution()

    def read_solution(self):
        """ read flux distribution/MFM from file """
        for solution in self.filenames:
            if self.method in ["fba", "moma"]:
                flux = MetabolicFlux()
                try:
                    flux.readFromFile(solution)
                except IOError as strerror:
                    print ("An error occurred while trying to read file %s:" %
                           os.path.basename(solution))
                    print(strerror)
                    exit()
                except SyntaxError as strerror:
                    print ("An error occurred parsing file %s:" %
                           os.path.basename(solution))
                    print(strerror)
                    exit()

                for rea in flux.fluxDict:
                    self.result[rea] = self.result[rea] + [[flux.fluxDict[rea]]] \
                        if rea in self.result else [[flux.fluxDict[rea]]]

            elif self.method == "fva":
                try:
                    flux, minmax = \
                        FvAnalyzer.parseSolutionFile(solution)
                except IOError as strerror:
                    print ("An error occurred while trying to read file %s:" %
                           os.path.basename(solution))
                    print(strerror)
                    exit()
                except SyntaxError as strerror:
                    print ("Error in FVA solution file %s:" %
                           os.path.basename(solution))
                    print(strerror)
                    exit()

                for rea in flux.fluxDict:
                    self.result[rea] = self.result[rea] + [[flux[rea]]+list(minmax[rea])] \
                        if rea in self.result else [[flux[rea]]+list(minmax[rea])]

            elif self.method == "mfm":
                solution = MfmAnalyzer.parseSolutionFile(solution)
                for met in solution:
                    self.result[met] = self.result[met] + [[solution[met]]] \
                        if met in self.result else [[solution[met]]]

            else:
                print ("Error: only FBA/MOMA, FVA, and MFM solutions can be read from file.\nIf "
                       "you want to calculate MVA, please provide two scenario files and a model!")
                exit()

    def is_sce(self):
        """ check if the given files are parameter files """
        for filename in self.filenames:
            try:
                # Parse file, get maxmin, name of objective function, and solver name
                pparser = ParamParser()
                pparser.parse(filename)
            except Exception:
                return False
        return True

    def do_analysis(self):
        """ perform FBA/FVA/MFM/MVA depending on given method """
        tolerance = 0.95

        for filename in self.filenames:
            model = MetabolicModel()
            rparser = ReactionParser()
            try:
                model.addReactionsFromFile(self.model, rparser)
            except IOError as strerror:
                print ("An error occurred while trying to read file %s:" %
                       os.path.basename(self.model))
                print(strerror)
                exit()
            except SyntaxError as strerror:
                print ("Error in reaction file %s:" %
                       os.path.basename(self.model))
                print(strerror)
                exit()
            # model = self.model
            model_messages = []
            pparser = ParamParser()
            # Parse file, get maxmin, name of objective function, and solver name
            maxmin, obj_str, solver, num_iter, lb_vec, ub_vec = pparser.parse(
                filename)
            fba_params = FbaParam(solver, maxmin, obj_str, num_iter)
            fba_params.setLinConstraints(pparser.lin_constraints)
            print(filename)
            # Set flux bounds in model
            model.setFiniteBounds(lb_vec, ub_vec, True, model_messages)

            solver = "default" if not solver else solver

            if self.method in ["fba", "moma"]:
                fba = FbAnalyzer(solver)
                _, solution, _ = fba.runOnModel(
                    model, fba_params, True, False, True)
                for rea in model.reactionDict:
                    self.result[rea] = self.result[rea] + [[solution[rea]]] \
                        if rea in self.result else [[solution[rea]]]

            elif self.method == "fva":
                fva = FvAnalyzer(solver)
                minmax, solution, _, _, _ = fva.runOnModel(
                    model, fba_params, tolerance)

                if not minmax:
                    print(
                        "FVA Warning: Model is infeasible or unbounded. Nothing to do.")
                    exit()
                for rea in model.reactionDict:
                    i = model.reactionDict[rea]
                    self.result[rea] = self.result[rea] + [[solution[i]]+list(minmax[i])] \
                        if rea in self.result else [[solution[i]]+list(minmax[i])]

            elif self.method == "mfm":
                mfm = MfmAnalyzer(solver)
                solution, _ = mfm.runOnModel(
                    model, fba_params, tolerance, minimize=True)
                for met in model.metaboliteDict:
                    i = model.metaboliteDict[met]
                    self.result[met] = self.result[met] + [[solution[i]]] \
                        if met in self.result else [[solution[i]]]

            elif self.method == "mva":
                fba = FbAnalyzer(solver)
                _, solution, _ = fba.runOnModel(
                    model, fba_params, True, False, True)

                mfm = MfmAnalyzer(solver)
                solution_min, _ = mfm.runOnModel(
                    model, fba_params, tolerance, minimize=True)
                solution_max, _ = mfm.runOnModel(
                    model, fba_params, tolerance, minimize=False)

                split_ratios = solution.computeAllSplitRatios(
                    model, 0., False, False)

                for met in model.metaboliteDict:
                    i = model.metaboliteDict[met]
                    flux = 0
                    for rea in split_ratios[met][0]:
                        flux += split_ratios[met][0][rea][1]

                    self.result[met] = self.result[met] + \
                        [[flux, solution_min[i], solution_max[i]]]\
                        if met in self.result else \
                        [[flux, solution_min[i], solution_max[i]]]

            else:
                print("No valid method was chosen.")
                exit()

    def set_plot_style(self, stylefile=None):
        """ read style for the scatterplot
            from configfile or take defaults """
        # first set defaults:
        config = configparser.SafeConfigParser(self.style)

        # try to read configfile, if exists
        if stylefile:
            try:
                open(stylefile)
            except IOError as error:
                print(error)
                print("Please enter a valid config file.")
                exit()
            config.read(stylefile)

        # parse style in a dictionary
        try:
            self.style = {
                "label_length": int(config.get("DEFAULT", "label_length")),
                "ids_in_names": config.get("DEFAULT", "ids_in_names") == "True",
                "unique_names": config.get("DEFAULT", "unique_names") == "True",
                "min_value": float(config.get("DEFAULT", "min_value")),
                "insignificant_color": config.get("DEFAULT", "insignificant_color"),
                "significant_color": config.get("DEFAULT", "significant_color"),
                "xlabel": config.get("DEFAULT", "xlabel"),
                "ylabel": config.get("DEFAULT", "ylabel"),
                "title": config.get("DEFAULT", "title"),
                "show_errorbar": config.get("DEFAULT", "show_errorbar"),
                "show_labels": config.get("DEFAULT", "show_labels"),
                "arrowcolor": config.get("DEFAULT", "arrowcolor"),
                "labelcolor": config.get("DEFAULT", "labelcolor"),
                "significance_threshold": float(config.get("DEFAULT", "significance_threshold")),
                "fontsize": float(config.get("DEFAULT", "fontsize")),
                "figsize": tuple(int(i) for i in config.get("DEFAULT", "figsize").split()),
                "plot_grid": config.get("DEFAULT", "plot_grid") == "True",
                "arrow_line_width": float(config.get("DEFAULT", "arrow_line_width")),
            }
        except Exception as error:
            print(error)
            print("Please check the configuration file")
            exit()

    def do_scatterplot(self, outputfile=None, obj_fun=None, stylefile=None):
        """ draw scatterplot """
        # read stylefile or get default style:
        self.set_plot_style(stylefile)
        style = self.style
        # calculate obj jun coefficient if obj fun is given
        if obj_fun:
            try:
                obj_coeff = old_div(self.result[obj_fun][0][0],
                                    self.result[obj_fun][1][0])
            except KeyError:
                print ("Error: Objective function %r could not be found in reactions dict!"
                       % obj_fun)
                exit()
        else:
            obj_coeff = 1

        # remove all entries that have only one result or are zero at some point
        self.result = {i: self.result[i] for i in self.result if len(self.result[i]) == 2
                       and self.result[i][0][0] > style["min_value"]
                       and self.result[i][1][0] > style["min_value"]}

        # create canvas for drawing the scatterplot
        plt.figure(figsize=style["figsize"])

        # parse result dict
        if self.method in ["fva", "mva"]:
            name, x, xmin, xmax, y, ymin, ymax = list(zip(*[
                ([rea] + self.result[rea][0] + self.result[rea][1]) for rea in self.result]))
        else:
            name, x, y = list(zip(*[([rea] + self.result[rea][0] + self.result[rea][1])
                                    for rea in self.result]))

        # calculate significance vector
        if self.method in ["fva", "mva"]:
            significance = [1 if abs_diff([xmin[i], xmax[i]],
                                          [ymin[i]*obj_coeff, ymax[i]*obj_coeff])
                            >= style["significance_threshold"]
                            else 0 for i in range(len(name))]
            # draw errorbar (only supported with fva and mva methods)
            if style["show_errorbar"] == "significant":
                plt.errorbar(x=[x[i] for i in range(len(x)) if significance[i] == 1],
                             y=[y[i]
                                 for i in range(len(y)) if significance[i] == 1],
                             xerr=[[xmin[i] for i in range(len(xmin)) if significance[i] == 1],
                                   [xmax[i] for i in range(len(xmax)) if significance[i] == 1]],
                             yerr=[[ymin[i] for i in range(len(ymin)) if significance[i] == 1],
                                   [ymax[i] for i in range(len(ymax)) if significance[i] == 1]],
                             ls='none', elinewidth=1, ecolor="black")
            elif style["show_errorbar"] == "all":
                plt.errorbar(x=x, y=y, xerr=[xmin, xmax], yerr=[ymin, ymax], ls='none',
                             elinewidth=1, ecolor="black")

        else:
            significance = [1 if old_div(abs(x[i]-y[i]*obj_coeff), (old_div((x[i]+y[i]*obj_coeff), 2)))
                            >= style["significance_threshold"] else 0 for i in range(len(name))]

        # add bisecting line
        plt.plot([0, 1000], [0, 1000], ls="--", c=".3")
        # format axis
        plt.ylim(min(y)-min(y)*0.1, (max(y)+max(y)*10))
        plt.xlim(min(x)-min(x)*0.1, (max(x)+max(x)*10))
        plt.yscale("log")
        plt.xscale("log")

        # add title and axis labels
        plt.xlabel(style["xlabel"])
        plt.ylabel(style["ylabel"])
        plt.title(style["title"])

        # add grid to plot (default is False)
        plt.grid(style["plot_grid"])

        # add scatter points
        plt.scatter(x, y, c=[style["significant_color"] if i == 1
                             else style["insignificant_color"] for i in significance])

        # add annotations
        name_abbr = abbreviate_all_in_list(name, max_len=style["label_length"],
                                           ids=style["ids_in_names"], unique=style["unique_names"])

        # use adjustText module if installed:
        if ADJ_TEXT_AV:
            texts = []
            for x_i, y_i, label, sig in zip(x, y, name_abbr, significance):
                if style["show_labels"] == "significant" and sig == 0:
                    continue
                elif style["show_labels"] not in ["all", "significant"]:
                    break
                texts.append(plt.text(x_i, y_i, label))
            if texts:
                adjust_text(texts, color=style["labelcolor"], fontsize=style["fontsize"],
                            arrowprops=dict(arrowstyle="->", color=style["arrowcolor"],
                                            lw=style["arrow_line_width"]))

        # otherwise just do annotations (overlap!)
        else:
            for i, label in enumerate(name_abbr):
                if style["show_labels"] == "significant" and significance[i] == 0:
                    continue
                elif style["show_labels"] not in ["all", "significant"]:
                    break
                plt.annotate(label, fontsize=style["fontsize"], xy=(x[i], y[i]), xytext=(6, 6),
                             color=style["labelcolor"], textcoords='offset points')
        # show interactive plot
        if outputfile:
            try:
                plt.savefig(outputfile, bbox_inches='tight', dpi=300)
                print("Scatterplot was saved as %r." % outputfile)
            except IOError:
                print("Error: No such file or directory: %r" % outputfile)
            except ValueError:
                print ("Error: The format of %r is not supported. "
                       "Please choose one of the following:\n"
                       "eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff"
                       % outputfile)
        else:
            plt.show()


def main():
    """main function"""

    version = "Scatterplot\n" + COPYRIGHT_VERSION_STRING
    parser = argparse.ArgumentParser(version=version)
    # required:
    parser.add_argument("-x", "--xfile", dest="xfile", required=True,
                        help="The input FILE for the x axis. Must be either FBA/FVA/MFM or a "
                        "scenario file for analysis.", metavar="FILE")
    parser.add_argument("-y", "--yfile", dest="yfile", required=True,
                        help="The input FILE for the y axis. Important: Must be of same kind "
                        "as the file for the x axis! ", metavar="FILE")
    parser.add_argument("-m", "--method", required=True,
                        dest="method", help="method for calculation or analysis. Choose "
                        "one of the following: fba, fva, mfm, mva (only with sce files)",
                        metavar="[FBA/FVA/MFM/MVA]")
    # optional:
    parser.add_argument("-r", "--reactions", dest="reactionfile", required=False, default=None,
                        help="Optional: perform Analysis on the network given "
                        "by the reaction FILE", metavar="FILE")
    parser.add_argument("-c", "--config", dest="configfile", required=False, default=None,
                        help="Optional: style file (ini) to define "
                        "the look of the scatterplot", metavar="FILE")
    parser.add_argument("-o", "--output", dest="outputfile", required=False, default=None,
                        help="Optional: write scatterplot to FILE",
                        metavar="FILE")
    parser.add_argument("-b", "--biomass", required=False, default=None,
                        dest="obj_fun", help="Optional: Name of objective function "
                        "to exclude from significance calculation ")

    options = parser.parse_args()

    rea = options.reactionfile
    # if rea:
    #     model = MetabolicModel()
    #     rparser = ReactionParser()
    #     try:
    #         model.addReactionsFromFile(rea, rparser)
    #     except IOError, strerror:
    #         print ("An error occurred while trying to read file %s:" %
    #                os.path.basename(rea))
    #         print strerror
    #         exit()
    #     except SyntaxError, strerror:
    #         print ("Error in reaction file %s:" %
    #                os.path.basename(rea))
    #         print strerror
    #         exit()

    scatter = Scatterplot(method=options.method,
                          filenames=[options.xfile, options.yfile],
                          model=rea if rea else None)

    scatter.do_scatterplot(outputfile=options.outputfile,
                           obj_fun=options.obj_fun,
                           stylefile=options.configfile)


if __name__ == '__main__':
    main()
