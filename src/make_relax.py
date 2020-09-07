from task import * 


class BWOMnRelax: 
      def __init__(self, workdir, encut, struct_path, name="relax_bwmn"):   

          """ 
          Sets input for Mn-doped BWO optimization  
          
          **Args:   
          """           
          potcar_path = "../pseudos/BWO_Mn_POTCAR" 
          kgrid = [2, 2, 2]   
          input_param =  DefaultOptimizationParameters(encut)  
          relax_calc = SCFCalculation(workdir, pseudo_par=None, kgrid=kgrid, name="BWO_Mn_relax", encut=encut, input_parameters=input_param)  
          relax_calc.make_calculation(struct_path, potcar_path=potcar_path)  

workdir_ = "/global/scratch/nleclerc/soc_mae_predictions/relax/BWO_Mn_doped/ortho/Mn2_Bi_site/v0" 
encut_ = 800 
struct_path_ = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/structures/BWO/BWO_doped/unrelaxed/BWO_orhto_Mn_Bi.vasp"   
test_relax = BWOMnRelax(workdir_, encut_, struct_path_)  
#class BWOMnRelax: 
#      def __init__(self, workdir, encut, struct_path, name="relax_bwmn"):   
#
#      """                                                                                                                                      
#      Sets input for Mn-doped BWO optimization  
#      
#      **Args:    
# 
#      """            
#      potcar_path = " "  
#      kgrid = [2, 2, 2]   
#      input_param =  DefaultOptimizationParameters(encut)  
#      relax_calc = SCFCalculation(workdir, pseudo_par=None, kgrid=kgrid, name="BWO_Mn_relax", encut=800, input_parameter=input_param)  
#      relax_calc.make_calculation(struct_path, potcar_path=potcar_path)  
