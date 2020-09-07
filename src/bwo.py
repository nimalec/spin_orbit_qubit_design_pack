from example_complexoxide import *

#npoints = 10   
#nbands_Fe_Bi = 1494
#workdir_Fe_Bi = "/global/scratch/nleclerc/soc_mae_predictions/bwo/bwoFe_Bi_orhto_mae_sphere_v01"
#struct_path_Fe_Bi = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/relaxed/BWO_orhto_Fe_Bi.vasp"  
#bwofe_orhto_Fe_Bi_mae = BWOFe_MAEworkflow(workdir_Fe_Bi, npoints, nbands_Fe_Bi, struct_path_Fe_Bi)
#bwofe_orhto_Fe_Bi_mae.make_calculations()

#npoints = 10 
#nbands_Co3_W = 1512 
#workdir_Co3_W = "/global/scratch/nleclerc/soc_mae_predictions/bwo/bwoCo3_W_ortho_sphere_v01"
#struct_path_Co3_W = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/relaxed/bwo_co3_w_pbesol_contcar.vasp" 
#
#bwo_orhto_Co3_W_mae = BWOCo_MAEworkflow(workdir_Co3_W, npoints, nbands_Co3_W, struct_path_Co3_W, nelect=1248)
#bwo_orhto_Co3_W_mae.make_calculations()    




#nbands_Mn4_W = 1512  
#workdir_Mn4_W = "/global/scratch/nleclerc/soc_mae_predictions/bwo/bwoMn4_W_ortho_sphere_v02"
#struct_path_Mn4_W = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/relaxed/bwo_mn4_w_pbesol_contcar.vasp"  
#
#bwo_orhto_Mn4_W_mae = BWOMn_MAEworkflow(workdir_Mn4_W, npoints, nbands_Mn4_W, struct_path_Mn4_W, nelect=1245)
#bwo_orhto_Mn4_W_mae.make_calculations()  




#nbands_Mn4_Bi = 1494 
#workdir_Mn4_Bi = "/global/scratch/nleclerc/soc_mae_predictions/bwo/bwoMn4_Bi_ortho_sphere_v01"
#struct_path_Mn4_Bi = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/relaxed/bwo_mn4_bi_pbesol_contcar.vasp" 
#
#bwo_orhto_Mn4_Bi_mae = BWOMn_MAEworkflow(workdir_Mn4_Bi, npoints, nbands_Mn4_Bi, struct_path_Mn4_Bi, nelect=1239)
#bwo_orhto_Mn4_Bi_mae.make_calculations()   

npoints = 200  
nbands_Fe3_bi = 1494 
workdir_Fe3_Bi = "/global/scratch/nleclerc/soc_mae_predictions/bwo/bwoFe3_Bi_pbca_sphere_v01"
struct_path_Fe3_Bi = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/relaxed/bwo_fe_bi_pbca_pbes.vasp" 
bwo_orhto_Fe3_Bi_mae = BWOFe_MAEworkflow(workdir_Fe3_Bi, npoints, nbands_Fe3_bi, struct_path_Fe3_Bi)
bwo_orhto_Fe3_Bi_mae.make_calculations()   
