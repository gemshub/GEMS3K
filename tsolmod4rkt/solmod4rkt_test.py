from solmod4rkt import *

task_data_file_name = "Thermo-time-in/series1-dat.lst"
task = SolModFactory(task_data_file_name)

print("Task: ", task_data_file_name)
print("PhSolNumber: ", task.Get_SolPhasesNumber())
print("PhSolNames: ", task.Get_SolPhasesNames())

# Getting SolModEngine for a feldspar phase by name
phase = task.SolPhase("Plagioclase");
print("Phase: name ", phase.Get_SolPhaseName(),
      ", mixing/activity model type ", phase.Get_MixModelType(),
      ", mixing/activity model code ", phase.Get_MixModelCode(),
      ", number of endmembers ", phase.Get_SpeciesNumber())

# Setting phase composition
x2m = {'Albite': 0.187, 'Anorthite': 3.5e-09, 'Sanidine': 0.813}
phase.SetMoleFractions(x2m)
print(phase.GetMoleFractions())

# Calculating activity coefficients of end members
phase.SolModActivityCoeffs()

# Writing results to a text file
phase.to_text_file("solmod_act_coef.txt", True)

# Get activity coefficients and print them in dict style
lnGamma2m = phase.GetlnActivityCoeffs()
print("Calculated activity coefficients: \n", lnGamma2m)


map_ideal = phase.SolModIdealProps()
print("Ideal properties of mixing: \n", map_ideal)

map_excess = phase.SolModExcessProps()
print("Excess properties of mixing: \n", map_excess)


