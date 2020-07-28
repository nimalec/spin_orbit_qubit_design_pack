from example_complexoxide import *   

workdir = "/global/scratch/nleclerc/soc_mae_predictions/bwoFe_orhto_mae_sphere_v02"   
npoints = 200 
nbands = 1494
struct_path = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO_orhto_Fe_Bi.vasp"    
bwofe_orhto_Fe_Bi_mae = BWOFe_MAEworkflow(workdir, npoints, nbands, struct_path) 
bwofe_orhto_Fe_Bi_mae.make_calculations()  
