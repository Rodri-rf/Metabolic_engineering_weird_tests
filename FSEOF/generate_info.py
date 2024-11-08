import cobra
from src.models import tb18
# import solvers included in optlang
import cobra.util.solver

"""general .sbml files for the different strains we are using
model1 = tb18.ecoli_pMeOH_pThr()
model3 = tb18.ecoli_BW25113()
model4 = tb18.ecoli_iML1515()

# export each model to a .sbml file
cobra.io.write_sbml_model(model1, "pMeOH.sbml")
cobra.io.write_sbml_model(model3, "BW25113.sbml")
cobra.io.write_sbml_model(model4, "iML1515.sbml")


#print the biomass reaction for each model
print(model1.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M"))
print(model3.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M"))
print(model4.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M"))
"""
# source: https://cobrapy.readthedocs.io/en/latest/solvers.html
print(cobra.util.solver.get_solver_name())
# print the solvers included in optlang
print(cobra.util.solver.solvers)
# run FBA on pMeOH.sbml once for each solver
model1 = cobra.io.read_sbml_model("pMeOH.sbml")
solution = model1.optimize()
print(solution.fluxes)
# change the solver to cplex
model1.solver = "cplex"
solution = model1.optimize()
print(solution.fluxes)
# change the solver to gurobi
model1.solver = "gurobi"
solution = model1.optimize()
print(solution.fluxes)
# change the solver to glpk
model1.solver = "glpk"
solution = model1.optimize()
print(solution.fluxes)
# change the solver to scipy
model1.solver = "scipy"
solution = model1.optimize()
print(solution.fluxes)


