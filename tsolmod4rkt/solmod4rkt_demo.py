import os
from solmod4rkt import *

# Path to the GEMS3K file set
task_data_file_name = "test01/gems3k-files/series1-dat.lst"

# Test file exist
if not os.path.exists(task_data_file_name):
    print( "File does not exist: ", task_data_file_name)
    exit(1)


# Initialize SolModFactory from the GEMS3K file set
task = SolModFactory(task_data_file_name)

print("Task:", task_data_file_name)
print("T(K): {} P(bar): {} N(PhSolutions): {}".format(
      task.Get_Temperature(), task.Get_Pressure(), task.Get_SolPhasesNumber()))

print("ElementNames:       \n", task.Get_AllElementNames())
print("ElementMoleAmounts: \n", task.Get_ElementMoleAmounts())

print("PhasesNames:      \n", task.Get_AllPhasesNames())
print("PhaseMoleAmounts: \n", task.Get_PhaseMoleAmounts())

print("SpeciesNames:       \n", task.Get_AllSpeciesNames())
print("SpeciesMoleAmounts: \n", task.Get_SpeciesMoleAmounts())

for k in range(1, task.Get_SolPhasesNumber()):
    phase = task.Sol_Phase(k)
    print("\nPhase: '{}'; mixing/activity model type: '{}'; model code: '{}'; N endmembers: {}".format(
          phase.Get_SolPhaseName(), phase.Get_MixModelType(), phase.Get_MixModelCode(), phase.Get_SpeciesNumber()))
    print("SpeciesNames:     ", phase.Get_SpeciesNames())
    phase.SolModActivityCoeffs()
    print("MoleFractions:    ", phase.Get_MoleFractions())
    print("lnActivities:     ", phase.Get_lnActivities())
    print("lnActivityCoeffs: ", phase.Get_lnActivityCoeffs())
    print(phase)


#print(task)
