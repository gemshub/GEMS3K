from solmod4rkt import *

task = SolModFactory("Thermo-time-in/series1-dat.lst")

phase = task.solution_phase("Plagioclase");
wx = {'Albite': 0.186993363098213, 'Anorthite': 3.45294711467247e-09, 'Sanidine': 0.81300663344884}
phase.SetMoleFractionsWx(wx, 0.)
print(phase)

phase.SolModActCoeff()
ln_gamma = phase.GetlnGamma()
print("\nln_gamma: ", ln_gamma)

map_ideal = phase.SolModIdealProp()
print("\nIdealProp: ", map_ideal)

map_excess = phase.SolModExcessProp()
print("\nExcessPropp: ", map_excess)

phase.SolModActCoeff()
ln_gamma = phase.GetlnGamma()
print("\nln_gamma: ", ln_gamma)

task.UpdateThermodynamic(973.15, 1000)
phase.SolModActCoeff()
ln_gamma = phase.GetlnGamma()
print("\nln_gamma: ", ln_gamma)

