from example_complexoxide import *

npoints = 200
#nbands_Fe_Bi = 1494
# workdir_Fe_Bi = "/global/scratch/nleclerc/soc_mae_predictions/bwoFe_orhto_mae_sphere_v02"
# struct_path_Fe_Bi = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO_orhto_Fe_Bi.vasp"
# bwofe_orhto_Fe_Bi_mae = BWOFe_MAEworkflow(workdir_Fe_Bi, npoints, nbands_Fe_Bi, struct_path_Fe_Bi)
# bwofe_orhto_Fe_Bi_mae.make_calculations()


nbands_Fe_W = 1728
workdir_Fe_W = "/global/scratch/nleclerc/soc_mae_predictions/bwoFe_W_orhto_mae_sphere_v0"
struct_path_Fe_W = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO_Fe_W_CONTCAR.vasp"
bwofe_orhto_Fe_W_mae = BWOFe_MAEworkflow(workdir_Fe_W, npoints, nbands_Fe_W, struct_path_Fe_W, nelect=1247) 
bwofe_orhto_Fe_W_mae.make_calculations()
