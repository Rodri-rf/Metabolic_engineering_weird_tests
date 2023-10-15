#!/usr/bin/env python

from multiprocessing import freeze_support
from cobra.io import read_sbml_model
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
from scipy.optimize import curve_fit
import sys
import argparse
from FSEOF_tests import FSEOF

# Parse arguments form the command line
parser = argparse.ArgumentParser(
    description="Identify gentetic targets for over-expression to increase the flux to a reaction of interest."
)

parser.add_argument("sbmlFile", type=str,
                    help="Path to the SBML model that should be scanned for over-expression targets. Needs to be a .xml file!")
parser.add_argument("biomassID", type=str, help="ID of the biomass reaction of the SBML file you provided")
parser.add_argument("reactionID", type=str,
                    help="ID of the reaction that will be optimized by genetic engineering. Over-expression targets will be identified for this reaction of the SBML model")

parser.add_argument("--iterations", type=int, action="store", default=30,
                    help="Number of iterations for the FSEOF algorithm. The default is 30. NOTE: This number should be decreased if flux variabillity is used.")
parser.add_argument("--useFVA", action="store_true",
                    help="Changes the method for finding over-expression from flux balance analysis to flux variabillity analysis. This will significantly increase the runtime from a few minutes to several hours!")
parser.add_argument("--constrainBiomass", action="store_true",
                    help="Constrains growth rate to 95%% of the theoretic maximum. This might improve the accuracy of the processed_results, but can also lead to mathmatical infeasible solutions that will cause an error.")
parser.add_argument("--changeBiomassConstrain", type=float, action="store", default=0.95,
                    help="If you would like to add an additional constrain to the growth rate, but not at 95%% of the theoretic maximum, use this option to specify at what percentage you want to set the constrain. NOTE: Percentages need to be passed as floats")
'''
#add argument to decide whether to stipulate a growth medium or not, it should be possible to pass a csv file with the medium
parser.add_argument("--medium", type=str, action="store", default=None,
                    help="Path to a .csv file containing the medium that should be used for the FSEOF algorithm. If no medium is provided, the default minimal medium will be used.")
'''
# add argument to decide whether to stipulate a growth medium or not. If argument is added,
# the minimal medium will be used and the argument will be two numbers: the lower and upper bound for methanol flux in the medium.
parser.add_argument("--medium", type=float, action="store", nargs=2, default=None,
                    help="Lower and upper bound for methanol flux in the medium. If no medium is provided, the default minimal medium will be used.")
# add argument to decide whether to minimize the target reaction
parser.add_argument("--minimize", action="store_false",
                    help="Minimize the target reaction instead of maximizing it. Default is to maximize the target reaction.")

args = parser.parse_args()


def main():
    f = FSEOF(args.sbmlFile, args.biomassID, args.reactionID, args.minimize)
    f.find_targets(args.steps, useFVA=args.useFVA, constrainBiomass=args.constrainBiomass,
                   maxFluxCutoff=args.changeBiomassConstrain, maximize=args.minimize)
    f.addReactionData()

    # Save reuslts. Depending on which method was used for the FSEOF algorithm, different information are exported.
    filename = args.reactionID

    if args.useFVA == True:
        f.targets.to_excel(
            "AmplificationTargets_{reaction}_FVA.xlsx".format(reaction=filename),
            columns=["q_slope", "l_sol", "q_slope_classifier", "l_sol_classifier", "reaction_class", "Reaction",
                     "Compartments"]
        )

    else:
        f.targets.to_excel(
            "AmplificationTargets_{reaction}.xlsx".format(reaction=filename),
            columns=["q_slope", "q_slope_Pearson correlation coefficient", "Reaction", "Compartments", "genes",
                     "target_reaction", "flux_change", "iterations", "strain"]
        )
        f.raw_data.to_excel("raw_data_FBA3.xlsx")


if __name__ == "__main__":
    freeze_support()  # can be neccessary for performing FVA to enable multiprocessing
    main()
