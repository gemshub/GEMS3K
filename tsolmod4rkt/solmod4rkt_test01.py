import math
from solmod4rkt import *

# Initialize SolModFactory from the GEMS3K file set
task_data_file_name = "Thermo-time-in/series1-dat.lst"
task = SolModFactory(task_data_file_name)

print("Task:", task_data_file_name)
print(" T(K): {} P(bar): {} N(PhSolutions): {}".format( task.Get_Temperature(), task.Get_Pressure(), task.Get_SolPhasesNumber()))
print("PhSolNames:", task.Get_SolPhasesNames())

# Getting SolModEngine for a feldspar phase 1 by name
phase1 = task.SolPhase("Alkali feldspar");
print("\nPhase1: name: '{}'; mixing/activity model type: '{}'; model code: '{}'; N endmembers: {}".format(
      phase1.Get_SolPhaseName(), phase1.Get_MixModelType(), phase1.Get_MixModelCode(), phase1.Get_SpeciesNumber()))

# Setting composition of the first feldspar phase (in mole fractions)
x1m = {'Albite': 0.20987, 'Anorthite': 1.7e-09, 'Sanidine': 0.79013 };
phase1.SetMoleFractions(x1m)

# Calculating activity coefficients of end members
phase1.SolModActivityCoeffs()

# Printing input phase composition and species activities
x_ph1 = phase1.GetMoleFractions()
a_ph1 = phase1.GetlnActivities()
for key in x_ph1:
    print("   '{}': x= {:.6g}; a= {:.6g}".format(key, x_ph1[key], math.exp(a_ph1[key])))

# Writing results to a text file
phase1.to_text_file("solmod_act_coef.txt", True)

# Get activity coefficients and print them
lnGamma1v = phase1.GetlnActivityCoeffs()
print("Calculated activity coefficients of endmembers:")
for key, value in lnGamma1v.items():
    print("   '{}': ln(gamma)= {:.6g}; gamma= {:.6g}".format(key, value, math.exp(value)))


# Getting SolModEngine for a feldspar phase 2 by index
phase2 = task.Sol_Phase(2);
print("\nPhase2: name: '{}'; mixing/activity model type: '{}'; model code: '{}'; N endmembers: {}".format(
      phase2.Get_SolPhaseName(), phase2.Get_MixModelType(), phase2.Get_MixModelCode(), phase2.Get_SpeciesNumber()))

# Setting composition of the second feldspar phase (in mole fractions)
x2m = {'Albite': 0.94371, 'Anorthite': 1.12e-07, 'Sanidine': 0.05629}
phase2.SetMoleFractions(x2m)

# Calculating activity coefficients of end members
phase2.SolModActivityCoeffs()

# Printing input phase 2 composition in dict style
print("  ", phase2.GetMoleFractions())

# Printing output activities
print("Calculated activities of endmembers: ")
for key, value in phase2.GetlnActivities().items():
    print("   '{}': a= {:.6g}".format(key, math.exp(value)))

# Writing results to a text file
phase2.to_text_file("solmod_act_coef.txt", True)

# Get activity coefficients and print them
lnGamma2v = phase2.GetlnActivityCoeffs()
print("Calculated activity coefficients of endmembers:")
for key, value in lnGamma2v.items():
    print("   '{}': ln(gamma)= {:.6g}; gamma= {:.6g}".format(key, value, math.exp(value)))


map_ideal = phase2.SolModIdealProps()
print("\nIdeal properties of mixing in phase2:\n", map_ideal)

map_excess = phase2.SolModExcessProps()
print("\nExcess properties of mixing in phase2:\n", map_excess)


