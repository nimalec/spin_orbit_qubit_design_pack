import numpy as np 
import os 
from os import path 

def get_total_energy_h(path):
   energ_list = []
   fl_nm = path+"/"+"OUTCAR" 
   isfile = os.path.isfile(fl_nm)
   if isfile == False:
        print("OUTCAR file not present! try to re-run the calculation.")
        pass
   else:
     with open(fl_nm, 'r') as f:
       for line in f.readlines():
         if 'TOTEN' in line:
           energ_list.append(line) 
   if len(energ_list) == 0: 
      return None  
   else:  
     tot_energ = energ_list[len(energ_list)-1] 
     return float(tot_energ[30:45])

def get_saxis(path):
    spin_lines = []  
    fl_nm = path+"/"+"INCAR" 
    isfile = os.path.isfile(fl_nm)     
    with open(fl_nm, 'r') as f:
        for line in f.readlines():
            if 'SAXIS' in line:
                spin = line 
                break 
    return spin[9:len(line)-1]      

def get_energy_list(workdir, dirs): 
    energ = []  
    line = workdir+"/"+"scf_ncl_" 
    for d in dirs:
        d_line = line + str(d) 
        toten = get_total_energy_h(d_line)
        energ.append(toten) 
    return energ  

def get_spin_list(workdir, dirs): 
    spins = []  
    line = workdir+"/"+"scf_ncl_" 
    for d in dirs:
        d_line = line + str(d) 
        s_line  = get_saxis(d_line)
        spins.append(s_line)  
    return spins 


def get_mae_list(workdir, dirs):  
    energ = get_energy_list(workdir, dirs)
    temp = []
    mae = []
    for i in energ:
        if i == None: 
           continue 
        else: 
           temp.append(i) 
    min_energ = min(temp) 
    for i in energ: 
        if i == None:  
           mae.append(None) 
        else: 
          mae.append(round((i-min_energ)*(10**(6)),4))   
    return mae 


def mae_outfile_write(workdir, dirs):
    flnm = workdir+"/""mae_output.txt"
    spins = get_spin_list(workdir, dirs) 
    energ = get_energy_list(workdir, dirs)   
    mae =  get_mae_list(workdir, dirs)   
    
    f = open(flnm, "w") 
    f.write("ID  SAXIS      TOTAL ENERGY (eV)     MAE (ueV)" + "\n\n") 
    for i in range(len(dirs)-1): 
        line = str(dirs[i])+"    "+spins[i]+"    "+str(energ[i])+"    "+str(mae[i])+"    "+"\n" 
        f.write(line) 
    f.close()             

