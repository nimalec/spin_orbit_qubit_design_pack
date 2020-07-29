import numpy as np
import os
import shutil
from structure import *

__author__ = 'Nima Leclerc'
__email__ = 'nleclerc@lbl.gov'

def make_runscript_h(work_dir, input_settings):
    """
    Generates runscript provided run settings for VASP calculation.

    **Args:

    workdir (str):
    input_settings (InputParameters):
    """
    run_settings = input_settings._parallel_settings
    fl_nm = run_settings["flnm"]
    if os.path.exists(work_dir+fl_nm) is True:
        os.remove(work_dir+fl_nm)

    f=open(fl_nm, "w")
    f.write("#!/bin/bash" + "\n\n")
    f.write("#SBATCH --job-name="+run_settings["job_name"]+"\n")
    f.write("#SBATCH --partition="+run_settings["partition"]+"\n")
    f.write("#SBATCH --account="+run_settings["machine"]+"\n")
    f.write("#SBATCH --qos=normal \n")
    f.write("#SBATCH --nodes="+str(run_settings["nodes"])+"\n")
    f.write("#SBATCH --ntasks-per-node=" + str(run_settings["ppn"])+"\n")
    f.write("#SBATCH --time="+run_settings["max_time"]+"\n\n")
    f.write("module unload intel/2016.4.072\n")
    f.write("module load intel/2018.5.274.par\n")
    f.write("module load vasp_intelmpi/5.4.4.16052018\n\n")
    f.write("EXE="+"'"+run_settings["exec"]+"'"+"\n\n")
    f.write("time mpirun $EXE\n\n")
    f.write("exit 0\n\n")
    f.close()
    shutil.move(fl_nm, work_dir)
    if os.path.exists("__pycache__") is True:
       os.system("rm -r __pycache__")

def make_incar_h(work_dir, input_settings,name="system"):
    """
    Generates VASP INCAR file for VASP calculation.

    **Args:

    workdir (str):
    input_settings (InputParameters):
    """

    fl_nm = "INCAR"
    if os.path.exists(work_dir+fl_nm) is True:
        os.remove(work_dir+fl_nm)
    else:
        pass
    f=open(fl_nm, "w")

    f.write("SYSTEM=   "+name+"\n")
    f.write("start parameters"+"\n")

    for key in input_settings._start_settings:
        if input_settings._start_settings[key]:
            f.write(key+"=   "+str(input_settings._start_settings[key])+"\n")
        else:
            pass

    f.write("\n")
    f.write("parallel settings"+"\n")
    if input_settings._parallel_settings["NCORE"] == None:
       pass  
    else:  
       f.write("NCORE"+"=   "+str(input_settings._parallel_settings["NCORE"])+"\n")
    if input_settings._parallel_settings["KPAR"] == None: 
       pass 
    else: 
       f.write("KPAR"+"=   "+str(input_settings._parallel_settings["KPAR"])+"\n\n")
    f.write("\n")
    f.write("electronic"+"\n")
    for key in input_settings._electronic_settings:
        if input_settings._electronic_settings[key]:
            f.write(key+"=   "+str(input_settings._electronic_settings[key])+"\n")
        else:
            pass
    f.write("\n")

    if input_settings._ionic_settings:
        f.write("ionic"+"\n")
        for key in input_settings._ionic_settings:
            if input_settings._ionic_settings[key]:
                f.write(key+"=   "+str(input_settings._ionic_settings[key])+"\n")
            else:
                pass
        f.write("\n")

    if input_settings._magnetic_settings:
        f.write("magnetic"+"\n")
        for key in input_settings._magnetic_settings:
            if input_settings._magnetic_settings[key]:
                if key == "SAXIS":
                    saxis = input_settings._magnetic_settings[key]
                    saxis_line = str(saxis[0])+" "+str(saxis[1])+" "+str(saxis[2])
                    f.write(key+"=   "+saxis_line+"\n")
                elif key == "MAGMOM":
                   line = " " 
                   magmom = input_settings._magnetic_settings[key] 
                   for i in magmom: 
                     line += str(i) + " "   
                   #magmom_line = str(magmom[0])+" "+str(magmom[1])+" "+str(magmom[2])  
                   f.write(key+"=   "+line+"\n")
                else:
                   f.write(key+"=   "+str(input_settings._magnetic_settings[key])+"\n")
            else:
                pass

    # if input_settings._hybrid_settings:
    #     f.write("hybrid"+"\n")
    #     for key in input_settings._hybrid_settings:
    #         if input_settings._hybrid_settings[key]:
    #             f.write(key+"="+str(input_settings._hybrid_settings[key])+"\n")
    #         else:
    #             pass
    #     f.write("\n")

    f.write("\n") 
    if input_settings._hubbard_settings:
        f.write("hubbard"+"\n")
        for key in input_settings._hubbard_settings:
            if input_settings._hubbard_settings[key]:
                if key == "LDAUL" or key == "LDAUJ" or key=="LDAUU":
                    line = ""
                    for i in input_settings._hubbard_settings[key]:
                        line += str(i)+ " "
                    f.write(key+"=   "+line+"\n")
                else:
                    f.write(key+"=   "+str(input_settings._hubbard_settings[key])+"\n")
            else:
                pass
        f.write("\n")

    if input_settings._misc_settings:
        f.write("misc"+"\n")
        for key in input_settings._misc_settings:
            if input_settings._misc_settings[key]:
                f.write(key+"=   "+str(input_settings._misc_settings[key])+"\n")
            else:
                pass
        f.write("\n")

    f.close()
    shutil.move(fl_nm, work_dir)
    if os.path.exists("__pycache__") is True:
       os.system("rm -r __pycache__")

def make_potcar_h(work_dir, pseudo_par):
    """
    Generates VASP POTCAR file for VASP calculation.

    **Args:

    workdir (str):
    pseudo_par (dict):
    """
    pseudo_dir = pseudo_par["directory"]
    pseudos = pseudo_par["flavor"]

    paths = []
    for pot in pseudos:
        pseudo_path = pseudo_dir + pot
        paths.append(pseudo_path)
    files = " ".join(paths)
    os.system("cat"+files+">> POTCAR")
    shutil.move("POTCAR", work_dir)

def make_poscar_h(work_dir, structure, number, species, name=None):
    """
    Generates VASP POSCAR file for VASP calculation.

    **Args:

    structure (Structure):
    workdir (str):
    """
    sites = structure._sites
    lattice = structure._lattice
    name = name or structure._name
    number = [str(num) for num in number]
    species_line = " ".join(species)
    number_line = " ".join(number)

    fl_nm = "POSCAR"
    if os.path.exists(work_dir+"/"+fl_nm) is True:
        os.remove(work_dir+"/"+fl_nm)

    f=open(fl_nm, "w+")
    f.write(name+"\n")
    f.write(str(1)+"\n")
    f.write(str(lattice[0][0])+"  "+str(lattice[0][1])+"  "+str(lattice[0][2])+"\n")
    f.write(str(lattice[1][0])+"  "+str(lattice[1][1])+"  "+str(lattice[1][2])+"\n")
    f.write(str(lattice[2][0])+"  "+str(lattice[2][1])+"  "+str(lattice[2][2])+"\n")
    f.write(species_line+"\n")
    f.write(number_line+"\n")
    f.write("Direct"+"\n")
    for site in sites:
        coord = site._coord
        coord_line = str(coord[0])+"  "+str(coord[1])+"  "+str(coord[2])+"\n"
        f.write(coord_line)

    f.close()
    shutil.move("POSCAR", work_dir)
    if os.path.exists("__pycache__") is True:
       os.system("rm -r __pycache__")

def make_kpoints_h(work_dir, kmesh, qshift=None):
    """
    Generates VASP POSCAR file for VASP calculation.

    **Args:
    workdir (str):
    kmesh (list):
    qmesh (list):
    """
    fl_nm = "KPOINTS"
    f=open(work_dir+"/"+fl_nm, "w+")
    f.write("Automatic mesh \n")
    f.write(str(0)+"\n") 
    f.write("Gamma"+"\n")
    f.write(str(kmesh[0])+"  "+str(kmesh[1])+"  "+str(kmesh[2])+"\n")
    if qshift:
       f.write(str(qmesh[0])+"  "+str(qmesh[1])+"  "+str(qmesh[2])+"\n")
    else:
        f.write(str(0)+"  "+str(0)+"  "+str(0)+"\n")
    f.close()
    if os.path.exists("__pycache__") is True:
       os.system("rm -r __pycache__")

#def make_mae_output_h() 
