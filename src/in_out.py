import numpy as np
import os
import shutil

def make_runscript_h(work_dir, input_settings):
    """
    Generates runscript provided run settings for VASP calculation.

    **Args:

    workdir (str):
    input_settings (InputParameters):
    """
    run_settings = input_settings._parallel_settings
    fl_nm = run_settings["flnm"]

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
    os.system("rm -r __pycache__")

def make_incar_h(work_dir, input_settings):
    """
    Generates VASP INCAR file for VASP calculation.

    **Args:

     workdir (str):
    input_settings (InputParameters):
    """

    fl_nm = "INCAR"
    f=open(fl_nm, "w")
    #f.write("SYSTEM="+  +"\n\n")
    f.write("start parameters"+"\n")

    for key in input_settings._start_settings:
        if input_settings._start_settings[key]:
            f.write(key+"=   "+str(input_settings._start_settings[key])+"\n")
        else:
            pass

    f.write("\n")
    f.write("parallel settings"+"\n")
    f.write("ncore"+"=   "+str(input_settings._parallel_settings["ncore"])+"\n")
    f.write("kpar"+"=   "+str(input_settings._parallel_settings["kpar"])+"\n\n")

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
                f.write(key+"="+str(input_settings._ionic_settings[key])+"\n")
            else:
                pass
        f.write("\n")


    if input_settings._magnetic_settings:
        f.write("magnetic"+"\n")
        for key in input_settings._magnetic_settings:
            if input_settings._magnetic_settings[key]:
                f.write(key+"="+str(input_settings._magnetic_settings[key])+"\n")
            else:
                pass
            f.write("\n")

    if input_settings._hybrid_settings:
        f.write("hybrid"+"\n")
        for key in input_settings._hybrid_settings:
            if input_settings._hybrid_settings[key]:
                f.write(key+"="+str(input_settings._hybrid_settings[key])+"\n")
            else:
                pass
        f.write("\n")

    if input_settings._hubbard_settings:
        f.write("hubbard"+"\n")
        for key in input_settings._hubbard_settings:
            if input_settings._hubbard_settings[key]:
                f.write(key+"="+str(input_settings._hubbard_settings[key])+"\n")
            else:
                pass
        f.write("\n")

    if input_settings._misc_settings:
        f.write("misc"+"\n")
        for key in input_settings._misc_settings:
            if input_settings._misc_settings[key]:
                f.write(key+"="+str(input_settings._misc_settings[key])+"\n")
            else:
                pass
        f.write("\n")

    f.close()
    shutil.move(fl_nm, work_dir)
    os.system("rm -r __pycache__")



def make_potcar_h(work_dir, pseudo_par):
    """
    Generates VASP POTCAR file for VASP calculation.

    **Args:

    workdir (str):
    pseudo_par (dict):
    """
    pseudo_dir = pseudo_par["dir"]
    pseudos = pseudo_par["flavor"]

    paths = []
    for pot in pseudos:
        pseudo_path = pseudo_dir + pot
        paths.append(pseudo_path)
    files = " ".join(paths)
    os.system("cat"+files+">> POTCAR")

# def make_poscar_h(work_dir, structure):
#     """
#     Generates VASP POSCAR file for VASP calculation.
#
#     **Args:
#
#     structure (Structure):
#     workdir (str):
#     """
