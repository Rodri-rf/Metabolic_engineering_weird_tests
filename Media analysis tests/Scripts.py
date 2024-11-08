import math
import cobra
import sys
# requires pip install XlsxWriter
import xlsxwriter
import os
import pandas
import pandas as pd
#import Callable from anotations so the IDE doesn't complain when writing the type of a function's parameters
from typing import Callable

sys.path.append("../")
from src import genes, media
from src.models import ecoli_BW25113
from cobra.medium import minimal_medium
import matplotlib.pyplot as plt


def setup_model() -> cobra.Model:
    """
    This function sets up the model for the simulation.
    It knocks out the frmA gene and adds the pUD9 plasmid.
    It then simulates the minimal media.
    :return:
    """
    model = ecoli_BW25113()
    genes.knockout_frmA(model)
    genes.add_pUD9(model)
    media.simulate_M9_minimal_media(model, aerobic=True, glucose=False, meoh=True, aas="all")
    # solution = model.optimize()
    return model


def print_exchange_reactions(model: cobra.Model) -> None:
    """
    This function prints the exchange reactions of the model
    :param model: cobrapy model
    model.medium is a dict where keys are names of fluxes of
    exchange reactions and the values are their upper bounds
    """

    for chemical, max_flux in model.medium.items():
        chem_name = (model.reactions.get_by_any(chemical))
        chem_name = chem_name[0].name
        print(f"id:{chemical}, max flux {max_flux}, reaction name: {chem_name}")

def generate_summary_for_exchange_fluxes(result_df: pd.DataFrame, model: cobra.Model, summary_df: pd.DataFrame) -> None:
    """
    This function generates a summary of the exchange fluxes.
    :param result_df: a pandas dataframe with the exchange fluxes, in a format such each row corresponds to a reaction
    and each row to the flux thtough that reaction, at different growth rates.
    :param model: the GSMM of our target strain
    :param summary_df: a pandas dataframe with the summary of the exchange fluxes

    """
    for curr_exchange_rect, data in result_df.iterrows():
        # index is a tuple and data is a pandas seriesº
        '''
        Series.describe(percentiles=None, include=None, exclude=None)
        Generate descriptive statistics, including those that summarize the central tendency, 
        dispersion and shape of a dataset’s distribution, excluding NaN values.
        '''
        print(f"Exchange reaction: {curr_exchange_rect}")
        stats = data.describe()  # stats is type pandas.core.series.Serie
        max_flux = stats["max"]
        min_flux = stats["min"]
        stdev_flux = stats["std"]

        # Now let us generate a list of growth rates for which the max and min occur
        # data[data == min_flux].index outputs another series; I just transform it into a list:
        min_rates = [x for x in data[data == min_flux].index]
        max_rates = [x for x in data[data == max_flux].index]

        print_mix_flux_data(data, max_flux, max_rates, min_flux, min_rates, stats, stdev_flux)
        # Now let us add this data to the summary dataframe
        summary_df = pd.concat([summary_df, pd.DataFrame({
            'Reaction': [model.reactions.get_by_id(str(curr_exchange_rect)).name],
            'Max Flux': [max_flux],
            'Min Flux': [min_flux],
            'Standard deviation': [stdev_flux],
            'Max Growth Rates': [', '.join(map(str, max_rates))],
            'Min Growth Rates': [', '.join(map(str, min_rates))]
        })], ignore_index=True)
    return summary_df


def simulate_minimal_media_robust_fluxes(model: cobra.Model, min_target: int, max_target: int)->None:
    '''
    The function minimal_medium by default obtains the medium with the lowest total import flux.
    This function needs two arguments: the model and the minimum growth rate (or other objective)
    the model has to achieve.

    This function iterates through different target growth rates
    min_target: minimum biomass flux. Given as float, multiplied by 10 for convenience
    max_target: minimum biomass flux. Given as float, multiplied by 10 for convenience
    '''
    min_target = math.ceil(min_target * 10)
    max_target = math.ceil(max_target * 10)
    # Run the FBA to minimize input fluxes and store result as a data frame
    result_df = minimal_medium(model, min_target / 10)
    result_df = result_df.to_frame(str(min_target / 10))  # turn to dataframe
    # runs all the simulations
    for i in range(min_target + 1, max_target + 1):
        min_growth_rate = i / 10
        result_df[str(min_growth_rate)] = minimal_medium(model, min_growth_rate)
    summary_df = pd.DataFrame(
        columns=['Reaction', 'Max Flux', 'Min Flux', 'Standard deviation', 'Max Growth Rates', 'Min Growth Rates'])
    result_df.to_csv("C:\\Users\\rodcs\\Desktop\\output3.csv")

    summary_df = generate_summary_for_exchange_fluxes(result_df, model, summary_df)
    summary_df.to_csv("C:\\Users\\rodcs\\Desktop\\summary_exchange_fluxes_FBA.csv")
        # Alternative: path in local machine C:\Users\rodcs\Desktop\results_minimal_media

def print_mix_flux_data(data, max_flux, max_rates, min_flux, min_rates, stats, stdev_flux):
    """
    Prints the data of the fluxes of the exchange reactions for different growth rates
    """
    print("MAX: {a} for growth rates {b}".format(a=max_flux, b=max_rates))
    print("MIN: {a} for growth rates {b}".format(a=min_flux, b=min_rates))
    print("Standard deviation: {a}".format(a=stdev_flux))
    print(data[data == (stats["max"])].index[0])

def plot_fluxes(result_df:pandas.DataFrame)->None:
    """
    Plots the fluxes of the exchange reactions for different growth rates
    :param result_df:
    :return: None
    """
    result_df.plot()
    plt.xticks(range(len(result_df.index)), result_df.index)  # Show all x-axis labels
    plt.tight_layout()  # Adjust spacing between labels
    plt.show()

def simulate_media_robust_fluxes(model: cobra.Model, min_growth_rate, max_gowth_rate,
                                algorithm: Callable[[cobra.Model], None] = None,
                                filename: str = "Media_analysis_results.xlsx") -> None:
    """
    model: a cobra model of the target bacteria strain
    min_growth_rate: minimal growth rate of the bacteria we want to start iterating with
    max_gowth_rate: maximal growth rate of the bacteria we want to iterate to
    algorithm: the algorithm we want to use to solve for the fluxes of the model. By default it is set to regular FBA
    return: None, writes a CSV file with the results
    Function that simulates the fluxes of the exchange reactions for different growth rates.
    This function iterates through different target growth rates and stores the results in a CSV file.
    It enfornces a growth rate by setting the biomass flux as the objective function and its upper bound and
    lower bound nearly equal to the target growth rate. The function then runs the algorithm and stores the fluxes
    of the exchange reactions in a dictionary. The dictionary is then converted into a dataframe and stored in a CSV file.
    """
    max_gowth_rate = math.ceil(max_gowth_rate * 10)
    min_growth_rate = math.ceil(min_growth_rate * 10)
    dic_fluxes = {}  # dictionary to store the fluxes of the exchange reactions

    """ First, let's get the raw data for the FBA and store it in a dict. The simulations are run
    for different "enforced" growth rates."""
    for curr_rate in range(min_growth_rate, max_gowth_rate + 1):
        print(f"current enforced flux is {curr_rate / 10}")
        model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').lower_bound = curr_rate / 10
        model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').upper_bound = (curr_rate + 1) / 10
        if(algorithm is None):
            solution = model.optimize() #The Model.optimize() function will return a Solution object
        else:
            solution = algorithm(model)

        for curr_ex_reaction in model.exchanges: #TODO should we also consider demand and sink reactions?
            #the keys are the reaction IDs, the values are lists of fluxes
            curr_reaction_ID = curr_ex_reaction.id
            curr_flux = solution.fluxes[str(curr_ex_reaction.id)]
            #if the reaction is not yet in the dictionary, add it with the default value of an empty list, then append
            #the current flux to the list
            dic_fluxes.setdefault(curr_reaction_ID, []).append(curr_flux)

    result_df_FBA = save_flux_data_to_dataframe(dic_fluxes, min_growth_rate, max_gowth_rate)

    #plot_fluxes(result_df_FBA)
    summary_df = pd.DataFrame(columns=['Reaction', 'Max Flux', 'Min Flux',
                                       'Standard deviation', 'Max Growth Rates', 'Min Growth Rates'])
    summary_df = generate_summary_for_exchange_fluxes(result_df_FBA, model, summary_df)
    # Construct the full file path using the 'os' module
    file_path = os.path.join(os.path.expanduser("~"), "Desktop", filename)
    # Write the two dataframes to two CSV files, located in file_path
    result_df_FBA.to_csv(file_path + "_raw.csv")
    summary_df.to_csv(file_path + "_summary.csv")

    """ writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
   #write both dataframes to the same file, but in different sheets

    result_df_FBA.to_excel(writer, sheet_name='Exchange_fluxes_FBA_raw', index=True)
    summary_df.to_excel(writer, sheet_name='Exchange_fluxes_FBA_summary', index=True)
    """

def save_flux_data_to_dataframe(flux_data: dict, min_growth_rate: int, max_growth_rate: int)->pandas.DataFrame:
    """
    Saves the flux data to a pandas DataFrame. The columns are the exchange reactions,
    and the rows are the fluxes at different growth rates.
    Specifically, the indexes are the names of the exchange reactions.
    :param flux_data:
    :return: dataframe with the flux data
    """
    column_headers = [x / 10 for x in range(min_growth_rate, max_growth_rate + 1)]
    print(flux_data)
    #turn dictionary into a dataframe
    result_df_FBA = pd.DataFrame.from_dict(flux_data, orient='index', columns=column_headers)
    result_df_FBA.index.name = 'Reaction'
    result_df_FBA.columns.name = 'Growth rate'

    return result_df_FBA


if __name__ == "__main__":

    # our e.coli
    model = setup_model()
    simulate_media_robust_fluxes(model, 0.1, 1.5, cobra.flux_analysis.pfba, "Media_analysis_results_pFBA3.xlsx")
    simulate_media_robust_fluxes(model, 0.1, 1.5, filename="Media_analysis_results_FBA3.xlsx")


    #IMPORTANT: the attribute "objective value" is defined differently for FBA and pFBA, SO BE CAREFUL.
    #That's why I used the "fluxes" attribute to get the biomass flux.
    solution = model.optimize()
    print(f"flux under FBA, all amino acids, for 'vanilla'engineered e.coli: {solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']}")
    pfba_solution = cobra.flux_analysis.pfba(model)
    print(f"flux under p_FBA, all amino acids, for 'vanilla'engineered e.coli on p_FBA: {pfba_solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']}")
    #Now, let's simulate the fluxes for different media conditions: namely no amino acids:
    media.simulate_M9_minimal_media(model, aerobic=True, glucose=False, meoh=True, aas="none")
    solution = model.optimize()
    print(f"flux for 'vanilla'engineered e.coli, no amino acids: {solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']}")
    pfba_solution = cobra.flux_analysis.pfba(model)
    print(f"flux for 'vanilla'engineered e.coli, no amino acids, on p_FBA: {pfba_solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']}")

    print("---------------------------------------------------------------------------")


    # Wild type for comparison
    model = cobra.io.read_sbml_model('../data/iBWG_1329.xml')
    media.simulate_M9_minimal_media(model, aerobic=True, glucose=False, meoh=True, aas="all")
    solution = model.optimize()
    print(f"flux for wild type e.coli, , all amino acids: {solution.objective_value}")
    # Now, let's simulate the fluxes for different media conditions: namely no amino acids:
    media.simulate_M9_minimal_media(model, aerobic=True, glucose=False, meoh=True, aas="none")
    solution = model.optimize()
    print(f"flux for wild type e.coli, no amino acids: {solution.objective_value}")

