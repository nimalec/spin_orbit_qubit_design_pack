from example_complexoxide import *

npoints = 200
nbands_Fe_Ti = 1536
workdir_Fe_Ti = "/global/scratch/nleclerc/soc_mae_predictions/ptoFe_W_tetrag_mae_sphere_v0"
struct_path_Fe_Ti = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/PTO/POSCAR_FePTO_tetrag"
ptofe_tetrag_Fe_Ti_mae = PTOFe_MAEworkflow(workdir_Fe_Ti, npoints, nbands_Fe_Ti, struct_path_Fe_Ti, nelect=1185)
ptofe_tetrag_Fe_Ti_mae.make_calculations()
