import pandas
from cobra.io import read_sbml_model
import cobra
import pandas as pd
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
import scipy
import sys
sys.path.append("../../")
from src import media
# import solvers included in optlang
import cobra.util.solver
from scripts import Sanity_check

class FSEOF():

    def __init__(self, path, biomassID, objectiveID, maximize = True):

        # Declare class variables that require user input
        self.model = read_sbml_model(path)
        self.objectiveID = objectiveID
        self.biomassID = biomassID

        # Declare other class variables
        self.initial_fluxes = None
        self.maxProductFlux = None
        self.minProductFlux = None
        self.initialProductFlux = None
        self.optimalGrowth = None
        self.solution = None
        self.targets = None
        self.sorted_targets = None
        self.raw_data = None

        # Calculate initial fluxes and maximal product flux
        self.calculate_intitial_fluxes()
        if maximize:
            self.calculate_maxmimal_product_flux()
        else:
            self.calculate_minimal_product_flux(maximize = maximize)

        # define e a media for the model:
        # media.simulate_M9_minimal_media(self.model, aerobic=True, glucose=False, meoh=True, aas="all")
        # define a meoh flux in mmol/gDW/h
        media.simulate_M9_minimal_media_MeOH_variablity(self.model, aerobic=True,
                                                        glucose=False, meoh_flux=10, aas="none")


    def calculate_intitial_fluxes(self):

        """
        This function calculates the initial fluxes of the model. It stores the optimal growth rate of the organism
        and the flux through the reaction of interest.
        """
        self.initial_fluxes = self.model.optimize()

        # store the initial flux through the reaction of interest
        self.initialProductFlux = self.initial_fluxes.fluxes[self.objectiveID]
        self.optimalGrowth = self.initial_fluxes.objective_value
        # sanity check
        Sanity_check.sanity_check(self.model)


    def calculate_maxmimal_product_flux(self):

        """
        Calculate & store the theoretical maximal flux through the reaction of interest
        """
        if self.maximize:
            self.model.objective = self.model.reactions.get_by_id(self.objectiveID)
            self.solution = self.model.optimize()
            self.model.objective = self.model.reactions.get_by_id(self.biomassID)
            self.maxProductFlux = self.solution.objective_value
        else:
            # if the goal is to minimize the flux through the reaction of interest, the maximal flux is the initial flux
            self.maxProductFlux = self.initialProductFlux

    def calculate_minimal_product_flux(self, maximize = False):

            """
            Calculate & store the theoretical minimal flux through the reaction of interest
            """

            if  maximize:
                # if the goal is to maximize the flux through the reaction of interest, the minimal flux is the initial flux
                self.minProductFlux = self.initialProductFlux
            else:
                self.model.objective = self.model.reactions.get_by_id(self.objectiveID)
                # change the setup to minimize the flux through the reaction of interest
                self.model.objective.direction = "min"
                self.solution = self.model.optimize()
                self.model.objective = self.model.reactions.get_by_id(self.biomassID)
                self.minProductFlux = self.solution.objective_value

    def find_targets(self, steps: int, maxFluxCutoff=0.95, useFVA=True, constrainBiomass=False, maximize = True):

        """
        DESCRIPTION:
        Find over-expression targets for the reaction of interest using the flux scanning based on enforced objective flux (FSEOF) algorithm.

        Either FBA or FVA can be used. FVA will significantly increase the runtime.

        PARAMETERS:
        iterations: int
            Number of iterations used for the FSEOF algorithm

        constrainBiomass: bool, default: False
            Add a constrain for the biomass reaction during the FSEOF algorithm. The constrain can further be modified with the maxFluxCutoff parameter.

        maxFluxCutoff: float, default: 0.95
            Percentage of the theoretical maximal biomass flux that should be enfored during the FSEOF algorithm, if the constrainBiomass parameter is set to True.

        RETURNS:
        None
        """

        fluxes = pd.DataFrame()
        fluxCapacity = pd.DataFrame()
        enforcedFluxes = list()

        # define helper functions for linear regression
        def get_lin_regresion_info(y):
            # try replacing with linear regression from scipy
            # which is scipy.stats.linregress(x, y=None, alternative='two-sided')
            lin_regresion_result = scipy.stats.linregress(enforcedFluxes, y)
            # popt, _ = curve_fit(f, enforcedFluxes, y)
            return lin_regresion_result

        if constrainBiomass == True:
            # constrain growth rate to defined percentage of optimal growth rate
            growthConstraint = self.model.problem.Constraint(
                self.model.reactions.get_by_id(self.biomassID).flux_expression, lb=maxFluxCutoff * self.optimalGrowth
            )
            self.model.add_cons_vars(growthConstraint)

        if useFVA == True:
            pass
            # Here be dragons. The FVA code was just distracting, but I'll put it back later.

        elif useFVA == False:
            # perform FSEOF algorithm with FBA
            raw_FBA_results = pd.DataFrame()
            for i in range(steps):
                # enforced_flux_value is the desired flux for the reaction of interest
                if maximize == True:
                    enforced_flux_value = self.initialProductFlux + (i / steps) * (
                                self.maxProductFlux - self.initialProductFlux)
                    print("enforced flux value for {}: ".format(self.objectiveID), enforced_flux_value)
                    print("initial product flux: ", self.initialProductFlux)
                    print("max product flux: ", self.maxProductFlux)


                elif maximize == False:
                    # enforce a decreasing flux through the target reaction
                    enforced_flux_value = self.initialProductFlux - (i / steps) * (
                                self.initialProductFlux - self.minProductFlux)
                    print("enforced flux value for {}: ".format(self.objectiveID), enforced_flux_value)

                # add the constraint to the model
                reaction = self.model.reactions.get_by_id(self.objectiveID)
                # Create a constraint to enforce the desired flux value
                constraint = self.model.problem.Constraint(
                    reaction.flux_expression,
                    ub=enforced_flux_value)
                self.model.add_cons_vars(constraint)
                enforcedFluxes.append(enforced_flux_value)

                # print the status of the variabels
                print("initial product flux: ", self.initialProductFlux)
                print("min product flux: ", self.minProductFlux)
                # change solver
                # I've done
                solution = self.model.optimize()
                # enfored flux value is the computed flux through the reaction of interest
                print("enforced flux value after FBA: ", solution.fluxes[self.objectiveID])
                # add fluxes to fluxes dataframe. Each column represents one step of the FSEOF algorithm. The column
                # name is the enforced flux.
                fluxes["Average_{}".format(i)] = solution.fluxes
                raw_FBA_results["Average_{}".format(i)] = solution.fluxes

            fluxes["q_slope"] = fluxes.apply(lambda y: get_lin_regresion_info(y).slope, axis=1)
            # add the target reaction that FSEOF is trying to optimize to the fluxes dataframe
            fluxes["target_reaction"] = self.objectiveID
            # add the pearson correlation coefficient of the linear regression to the fluxes dataframe
            fluxes["q_slope_Pearson correlation coefficient"] = raw_FBA_results.apply(lambda y: (get_lin_regresion_info(y).rvalue), axis=1)
            fluxes["strain"] = self.model.id
            fluxes["iterations"] = steps
            # for each reaction also list its associated gene(s)
            fluxes["genes"] = fluxes.apply(lambda y: self.model.reactions.get_by_id(y.name).gene_reaction_rule, axis=1)
            print(self.model.medium)

            # calculate the change in flux through the target reaction from the first to the last step of the FSEOF algorithm
            fluxes["flux_change"] = raw_FBA_results.apply(lambda y: y["Average_{}".format(steps - 1)] - y["Average_0"], axis=1)
            # calculate the change in the flux through the biomass reaction from the first to the last step of the FSEOF algorithm
            self.targets = fluxes
            self.raw_data = raw_FBA_results

    def sort_targets(self, useFVA=False):
        if useFVA == True:
            # Here be dragons. The FVA code was just distracting, but I'll put it back later.
            pass
        self.targets.sort_values("q_slope", ascending=False, inplace=True)

    def addReactionData(self):
        metabolites = list()
        compartments = list()

        for reactionID in list(self.targets.index):
            metabolites.append(
                self.model.reactions.get_by_id(reactionID).build_reaction_string(use_metabolite_names=True))
            compartments.append(self.model.reactions.get_by_id(reactionID).get_compartments())

        self.targets["Reaction"] = metabolites
        self.targets["Compartments"] = compartments







