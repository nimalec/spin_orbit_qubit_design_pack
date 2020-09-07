from example_complexoxide import *

npoints = 160
nbands = 1428
workdir = "/global/scratch/nleclerc/soc_mae_predictions/pto/ptoFe_Ti_tetrag_sphere_v01"
struct_path = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/PTO/PTO_Fe/POSCAR_tetrag" 
pto_calc = PTOFe_MAEworkflow(workdir, npoints, nbands, struct_path,  nelect=1185) 
pto_calc.make_calculations()   
