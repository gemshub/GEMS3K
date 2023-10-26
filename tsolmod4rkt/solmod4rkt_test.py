from solmod4rkt import *

task = SolModFactory("Thermo-time-in/series1-dat.lst")

phase = task.SolPhase("Plagioclase");
wx = {'Albite': 0.186993363098213, 'Anorthite': 3.45294711467247e-09, 'Sanidine': 0.81300663344884}
phase.SetMoleFractions(wx)
print(phase)

phase.SolModActivityCoeffs()
ln_gamma = phase.GetlnActivityCoeffs()
print("\nln_gamma: ", ln_gamma)

map_ideal = phase.SolModIdealProps()
print("\nSolModIdealProps: ", map_ideal)

map_excess = phase.SolModExcessProps()
print("\nSolModExcessProps: ", map_excess)

phase.SolModActivityCoeffs()
ln_gamma = phase.GetlnActivityCoeffs()
print("\nln_gamma: ", ln_gamma)

task.UpdateThermoData(973.15, 1000)
phase.SolModActivityCoeffs()
ln_gamma = phase.GetlnActivityCoeffs()
print("\nln_gamma: ", ln_gamma)

