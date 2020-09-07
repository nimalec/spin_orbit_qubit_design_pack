from spin_constrain import * 

npoints = 25 

workdir = "/global/scratch/nleclerc/soc_mae_predictions/spin_multiplets/bwo_fe3_bi/spin_v1" 
struct_path = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/relaxed/BWO_orhto_Fe_Bi.vasp"  
nupdown_rng = [0.5,5.5] 
spin_calc = BWOFeSpinConstr(workdir, npoints, nupdown_rng,struct_path)     
spin_calc.make_calculations() 
 
