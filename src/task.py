import numpy as np
import os
from shutil import copyfile
from structure import *
from in_out import *
import math 
__author__ = 'Nima Leclerc'
__email__ = 'nleclerc@lbl.gov'


class InputParameters:
        def __init__(self, name=None, start_settings=None,
        parallel_settings=None, electronic_settings=None, magnetic_settings=None,
        ionic_settings=None, hubbard_settings=None, hybrid_settings=None, misc_settings=None):

            """
            Sets standard input parameters for a VASP calculation.

            **Args:

            name (str): name of input paramter settings [default = "input_param"]
            start_settings (dict): start settings for calculation [default = {"nwrite": 2, "istart": 1, "iniwav": 1,
             "icharg": None, "nelect": None, "lorbit": 11,
              "nedos": 1000, "loptics": False, "lelf": None, "lvhar": None, "rwigs": None, "lvtof": None}]
            parallel_settings (dict): parallization settings for calculation [default = "flnm": "scf_run.sh", "job_name": "scf_std", "machine": "nano" ,
             "partition": "etna", "nodes": 4,"ppn": 24,
              "max_time": 24.0, "ncore": 8, "kpar": 2, "exec": "vasp_std"}]
            electronic_settings (dict): electronic structure settings for calculation [default = {"algo": "Normal", "encut": 800,
            "nelm": 200, "nelmin": 4, "ediff": 10E-5, "ismear": 1,
            "sigma": 0.2,"lasph": True, "lreal": "Auto", "addgrid": True, "maxmix": 100, "bmix": 1.5}]
            ionic_settings (dict): ionic settings for calculation, used for relaxations and structure optimization [default=None]
            magnetic_settings (dict): magnetic structure settings for calculation [default=None]
            hybrid_settings (dict): hybrid/hse settings for accurate band calculation [default=None]
            hubbard_settings (dict): hubbard calculation settings for localized d/f orbital predicitons [default=None]
            misc_settings_settings (dict): miscalaneous settings for calculation, can be any VASP setting [default=None]

            """

            self._name = name or "input_param"
            self._start_settings = start_settings or {"NWRITE": 2, "ISTART": 1, "INIWAV": 1,
             "ICHARG": None, "NELECT": None, "LORBIT": 11,
              "NEDOS": 1000, "LOPTICS": ".FALSE.","ISYM": -1 , "LELF": None, "LVHAR": None, "RWIGS": None, "LVTOF": None, "NBANDS": None, "LWAVE": None}
            self._parallel_settings = parallel_settings or {"flnm": "run_scf.sh", "job_name": "scf_std", "machine": "nano" ,
             "partition": "etna", "nodes": 4,"ppn": 24,
              "max_time": "24:00:00", "NCORE": 8, "KPAR": 2, "exec": "vasp_std"}
            self._electronic_settings = electronic_settings or  {"PREC":"Accurate" , "ALGO": "Normal", "ENCUT": 800,
            "NELM": None, "NELMIN": None, "GGA": "PS" ,"EDIFF": 10E-05, "ISMEAR": 0,
            "SIGMA": 0.2, "LASPH": ".TRUE.", "LREAL": "Auto", "ADDGRID": ".TRUE.", "MAXMIX": 100, "BMIX": 1.5}
            self._ionic_settings = ionic_settings
            self._magnetic_settings = magnetic_settings
            self._hybrid_settings = hybrid_settings
            self._hubbard_settings = hubbard_settings
            self._misc_settings = misc_settings

        def get_input_settings(self):

            """Getter function for Vasp Input paramters"""

            input_settings = {"name": name,
             "start": self._start_settings, "parallel": self._parallel_settings ,
            "electronic": self._electronic_settings, "magnetic": self._magnetic_settings,
            "hybrid": self._hybrid_settings, "hubbard": self._hubbard_settings, "misc_setting": self._misc_settings}
            return input_settings

        def update_start_settings(self, key, value):

            """
            Update a parameter in start settings.

            **Args:

            key (str): key in start settings
            value : value corresponding to updated key
            """

            if key in self._start_settings:
                self._start_settings[key] = value
            else:
                print("key does not exist!! keys include: {charge_option, prec, encut, nstep, epsilon, pseudo, n_elect.structure, smear, sigma, isym}")

        def update_parallel_settings(self, key, value):

            """
            Update a parameter in parallel settings.

            **Args:

            key (str): key in parallel settings
            value : value corresponding to updated key
            """

            if key in self._parallel_settings:
                self._parallel_settings[key] = value
            else:
                print("key does not exist!! keys include: {flnm , job_name , machine, partition, nodes  ,ppn, max_time , ncore,  kpar}")

        def update_electronic_settings(self, key, value):
            """
            Update a parameter in electronic settings.

            **Args:

            key (str): key in electronic settings
            value : value corresponding to updated key
            """

            if key in self._electronic_settings:
                self._electronic_settings[key] = value
            else:
                print("key does not exist!! keys include: {prec_level, algo, encut , nelm,nelmin, ediff, sigma, lasph, lreal, addgrid, bmaxmix, bmix}")

        def update_ionic_settings(self, key, value):
            """
            Update a parameter in ionic settings.

            **Args:

            key (str): key in ionic settings
            value : value corresponding to updated key
            """
            if self._ionic_settings:
              if key in self._ionic_settings:
                self._ionic_settings[key] = value
              else:
                print("key does not exist!! keys include: {ediff ,nsw, ibrion ,isif, isym, nblock,  kblock}")
            else:
              print("magnetic settings not present!")

        def update_magnetic_settings(self, key, value):

            """
            Update a parameter in magnetic settings.

            **Args:

            key (str): key in magnetic settings
            value : value corresponding to updated key
            """

            if self._magnetic_settings:
              if key in self._magnetic_settings:
                self._magnetic_settings[key] = value
              else:
                print("key does not exist!! keys include: {ispin, magmom, nupdown, saxis, lsorbit,noncollinear}")
            else:
              print("magnetic settings not present!")


        def update_hubbard_settings(self, key, value):

            """
            Update a parameter in Hubbard settings.

            **Args:

            key (str): key in Hubbard settings
            value : value corresponding to updated key
            """

            if self._hubbard_settings:
              if key in self._hubbard_settings:
                  self._hubbard_settings[key] = value
              else:
                  print("key does not exist!! keys include: {ldau, ldatype, ldaul, dlauu, ldauj, lmaxmix}")
            else:
               print("hybrid settings not present!")


class DefaultOptimizationParameters(InputParameters):
        def __init__(self, encut, name="relax_settings"):
            """
            Sets default input parameters for optimization

            **Args:

            encut (float): planewave energy cutoff for calculation
            name (str): name for relaxation setting [default="relax_settings"]

            """

            ionic = {"EDIFF": 1E-17, "NSW": 20, "IBRION": 2,"ISIF": 2, "ISYM": -1, "NBLOCK": 1,  "KBLOCK": 20}
            InputParameters.__init__(self, ionic_settings=ionic, name=name)
            self.update_electronic_settings("ENCUT", encut)



class DefaultSCFParameters(InputParameters):
         def __init__(self, encut, name="scf_settings"):
             """
             Sets default input parameters for scf ground state energy calculation

             **Args:
               encut (float): planewave energy cutoff for calculation
               name (str): name for scf setting [default="scf_settings"]

             """
             InputParameters.__init__(self, name=name)
             self.update_electronic_settings("ENCUT", encut)

class DefaultSCFUParameters(InputParameters):
         def __init__(self, encut, ldaul, Uparam, Jparam, name="DFTU_settings"):
             """
             Sets default input parameters for scf ground state energy calculation with +U correction

             encut (float): planewave energy cutoff for calculation
             ldaul (list): list of  orbital types for each species
             Uparam (list): list of U parameters for each species
             Jparam (list): list of J paramters for each species
             name (str):  name for scf+U setting [default="DFTU_settings"]

             """

             dftu_settings = {"LDAU": ".TRUE." , "LDAUU": Uparam, "LDATYPE": 2, "LADAUL": ldaul, "LDAUJ": Jparam , "LMAXMIX": 4}
             InputParameters.__init__(self, name=name, hubbard_settings=dftu_settings)
             self.update_electronic_settings("ENCUT", encut)


class DefaultMagCLParameters(InputParameters):
         def __init__(self, encut, magmom, ldaul, Uparam, Jparam, nupdown=None ,name="DFTCL_settings"):
             """
             Sets default input parameters for scf spin collinear calculation

             encut (flt): planewave energy cutoff for calculation
             magmom (list): list of magnetic moments for each species
             ldaul (list): list of  orbital types for each species
             Uparam (list): list of U parameters for each species
             Jparam (list): list of J paramters for each species
             nupdown (flt): difference between up and down spins   
             name (str):  name for magnetic noncolinear calculation setting [default="DFTCL_settings"]

             """

             cl_settings =  {"ISPIN": 2, "MAGMOM": magmom, "SAXIS": None, "LSORBIT": None, "LNONCOLLINEAR": None, "NUPDOWN": nupdown}
             dftu_settings = {"LDAU": ".TRUE.", "LDAUU": Uparam, "LDATYPE": 2, "LDAUL": ldaul, "LDAUJ": Jparam , "LMAXMIMX": 4}
             InputParameters.__init__(self, name=name, magnetic_settings=cl_settings, hubbard_settings=dftu_settings)
             self.update_electronic_settings("encut", encut)

class DefaultMagNCLParameters(InputParameters):
         def __init__(self, encut, spinaxis, ldaul, Uparam, Jparam, nupdown=None, name='DFTCL_settings'):
             """
            Sets default input parameters for scf spin non-collinear calculation

             encut (flt): planewave energy cutoff for calculation
             spinaxis (ndarray): spinaxis  for calculation
             ldaul (list): list of  orbital types for each species
             Uparam (list): list of U parameters for each species
             Jparam (list): list of J paramters for each species
             nupdown (flt): differnce between number of up and down spins 
             name (str):  name for magnetic noncolinear calculation setting [default="DFTNCL_settings"]
             """
             ncl_settings =  {"ISPIN": 2, "MAGMOM": None, "SAXIS": spinaxis, "LSORBIT": ".TRUE.", "LNONCOLLINEAR": ".TRUE.", "NUPDOWN":nupdown}
             dftu_settings = {"LDAU": ".TRUE.", "LDAUU": Uparam, "LDATYPE": 2, "LDAUL": ldaul, "LDAUJ": Jparam , "LMAXMIX": 4}
             InputParameters.__init__(self, name=name, magnetic_settings=ncl_settings, hubbard_settings=dftu_settings)
             self.update_electronic_settings("ENCUT", encut)

class SCFCalculation():
     def __init__(self, workdir, pseudo_par, kgrid, structure=None, name="scf_calc", encut=600, input_parameters=None):
         """
         Sets standard input parameters for a VASP calculation.

         **Args:

         workdir (str): workdirectory for scf calculation
         pseudo_par (dict): pseudopotential parameters {"directory":  , "flavor":[]}
         kgrid (list): kgrid for calculation
         structure (Structure): structure for scf calculation
         name (str): calculation name [default="scf_calc"]
         encut (float): planewave energy cutoff
         input_paramters (InputParameters): input paramters for scf calculation [defualt=DefaultSCFParameters(encut)]

         """

         self._name  = name
         self._workdir = workdir
         self._structure = structure
         self._kmesh = kgrid
         self._pseudo_par = pseudo_par
         self._input_settings = input_parameters or DefaultSCFParameters(encut=encut)
         self._struct_path = None
         self._potcar_path = None
         self._kpoint_path = None
         self._struct_path = None
         self._run_status = "unstarted"
         self._jobid = None
         self._cputime = None
         self._tot_energy = None
         self._fermi = None
# add Gamma line between 0 and k points

     def make_calculation(self, struct_path=None, run_script_path=None, k_points_path=None, potcar_path=None):
         """
         Sets up VASP input files and directory

         **Args:

         struct_path (str): path (including POSCAR file) of availible POSCAR in external directory
         run_script_path (str): path (including runscript file) of availible runscript in external directory
         k_points_path (str): path (including kpoints file) of availible kpoints file in external directory

         """

         os.mkdir(self._workdir)
         print("Work Directory now in: " + self._workdir)
         make_incar_h(self._workdir, self._input_settings)

         if struct_path:
            copyfile(struct_path, self._workdir+"/"+"POSCAR")
            self._struct_path = struct_path
         else:
             #make_poscar_h(self._workdir, self._structure, [4], ["Mn"])
             pass

         if run_script_path:
             copyfile(run_script_path, self._workdir+"/"+"run_scf.sh")
         else:
             make_runscript_h(self._workdir, self._input_settings)

         if k_points_path:
             copyfile(k_points_path, self._workdir+"/"+"KPOINTS")
             self._kpoint_path = k_points_path
         else:
             make_kpoints_h(self._workdir, self._kmesh)

         if potcar_path:
             copyfile(potcar_path, self._workdir+"/"+"POTCAR")
             self._potcar_path = potcar_path
         else:
             make_potcar_h(self._workdir, self._pseudo_par)

         if os.path.exists("__pycache__") is True:
            os.system("rm -r __pycache__")

     def run_calculation(self):
         os.system("sbatch"+" "+self._input_settings._parallel_settings["flnm"])
         self._run_status = "submitted"

     def get_run_time(self):
        fl_nm = self._workdir + 'OUTCAR'
        isfile = os.path.isfile(fl_nm)
        if isfile == False:
            print("OUTCAR file not present! try to re-run the calculation.")
            pass
        else:
           with open(fl_nm, 'r') as f:
             for line in f.readlines():
               if 'Total CPU time used (sec):' in line:
                   time_str = line
        return float(time_str[49:57])

     # def update_run_status(self):
     #     current_status_ = self._run_status
     #     path_ = self._workdir
     #
     #     def check_slurm_file_h(path):  ## checks if slurm file is present
     #         slurm_files = []
     #         slurm_file = None
     #         slurm_files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and "slurm" in i]
     #         str_len = len(slurm_files)
     #
     #         if str_len == 0:
     #             pass
     #         else:
     #             slurm_file = slurm_files[str_len-1]
     #         return slurm_file
     #
     #     def check_outcar_file_h(path): ## checks if job is done, conditioned on there being outcar
     #         assert os.path.exists("OUTCAR"), "OUTCAR file not present in directory"
     #         isdone = False
     #         with open(fl_nm, 'r') as f:
     #             for line in f.readlines():
     #                 if "General timing" in line:
     #                     isdone = True
     #                     break
     #        return isdone
     #
     #     def is_submitted_h(path, current_status): gets job status
     #         if current_status == "submitted" or current_status == "running" or current_status == "finished":
     #             status = True
     #         else:
     #             status = False
     #        return status
     #
     #     def is_running_h(path, current_status):
     #
     #         if current_status == "running":
     #             job_status = True
     #         elif current_status == "finished" or current_status == "not_submitted":
     #             job_status = False
     #         else:
     #             if check_slurm_file_h(path) == None:
     #                 status = "submitted"
     #             else:
     #                 status = "running"
     #
     #     def fetch_run_id_h(path):
     #         flnm = check_slurm_file_h(path)
     #         assert flnm != None, "job not yet submitted, ID does not exist!"
     #         return str(flnm[6:14])
     #
     #  slurm_file_status = check_slurm_file_h(path_)
     #  outcar_status = check_outcar_file_h(path_)
     #  if self._run_status == "not_submited":
     #     pass
     #  elif self._run_status == "submited" and slurm_file_status == None:
     #       pass
     #  elif self._run_status == "submited" and slurm_file_status != None:
     #       self._run_status == "running"
     #       self._jobid = fetch_run_id_h(path_)
     #  elif self._run_status == "running" and check_outcar_file_h(path) == True:
     #       self._cputime = self.get_run_time()

#def get_total_energy(self):
#    energ_list = []
#    fl_nm = self._workdir + 'OUTCAR'
#    isfile = os.path.isfile(fl_nm)
#    if isfile == False:
#      print("OUTCAR file not present! try to re-run the calculation.")
#       pass
#    else:
#      with open(fl_nm, 'r') as f:
#        for line in f.readlines():
#          if 'TOTEN' in line:
#            energ_list.append(line)
#    tot_energ = energ_list[len(energ_list)-1]
#    return float(tot_energ[30:40])
#
#def get_fermi(self):
#    fl_nm = self._workdir + 'OUTCAR'
#    isfile = os.path.isfile(fl_nm)
#    if isfile == False:
#        print("OUTCAR file not present! try to re-run the calculation.")
#        pass
#    else:
#       with open(fl_nm, 'r') as f:
#         for line in f.readlines():
#           if 'E-fermi :' in line:
#               fermi_str = line
#   return float(fermi_str[12:18])

class SpinMultipletSeries:
     def __init__(self, workdir, npoints, nupdown_rng ,kgrid, nbands, nodes, ppn, ref_orient, ldaul, magmom, Uparam, Jparam, encut, potcar_path, struct_path, ismear, saxis=None, name ="scf", time="12:00:00",  sigma=0.2, nelect=None): 

         """
         Computes the energy dependence on total configuration
         """
         self._workdir = workdir
         self._npoints = npoints
         self._name = name
         self._structure_path = struct_path
         self._potcar_path =  potcar_path
         self._nupdownlst = np.linspace(nupdown_rng[0], nupdown_rng[1], npoints)

     
         def set_calcs_h(nupdown, saxis, workdir, nodes, ppn, time,ismear, sigma, kgrid, encut, magmom, ldaul, Uparam, Jparam, nelect, nbands=None):  
             settings  = []  
             for s in nupdown.tolist():  
                 if saxis: 
                    ncl_settings = DefaultMagNCLParameters(encut=encut, spinaxis=saxis, ldaul=ldaul, Uparam=Uparam, Jparam=Jparam, nupdown=s)
                    ncl_settings.update_parallel_settings("exec", "vasp_ncl") 
                    ncl_settings.update_electronic_settings("EDIFF", 1.0E-4)
                    settings.append(ncl_settings)                    
                 else:  
                    cl_settings = DefaultMagCLParameters(encut=encut, magmom=magmom, ldaul=ldaul, Uparam=Uparam, Jparam=Jparam, nupdown=s)  
                    cl_settings.update_electronic_settings("EDIFF", 1.0E-6)
                    settings.append(cl_settings)
             calcs = []  
             itr = 0 
             for s in settings: 
                 s.update_start_settings("NBANDS", nbands)
                 s.update_start_settings("LWAVE", ".FALSE.")
                 s.update_parallel_settings("flnm ", "run.sh")
                 s.update_parallel_settings("job_name", "scf"+str(itr))
                 s.update_parallel_settings("nodes", nodes)
                 s.update_parallel_settings("ppn", ppn)
                 s.update_parallel_settings("max_time", time)
                 s.update_parallel_settings("KPAR", None)
                 s.update_electronic_settings("ISMEAR", ismear)
                 s.update_electronic_settings("SIGMA", sigma)
                 s.update_electronic_settings("EDIFF", 1.0E-6)
                 s_dir = workdir+"/"+"scf_"+str(itr)
                 calc = SCFCalculation(s_dir, pseudo_par=None, kgrid=kgrid, name="scf_"+str(itr), input_parameters=s) 
                 itr += 1 
                 calcs.append(calc) 
             return calcs 
 
         self._calc_list = set_calcs_h(self._nupdownlst, saxis, workdir, nodes, ppn, time,ismear, sigma, kgrid, encut, magmom, ldaul, Uparam, Jparam,  nelect, nbands) 
 
     def make_calculations(self):
         os.mkdir(self._workdir)
         for calc in self._calc_list: calc.make_calculation(struct_path=self._structure_path, potcar_path=self._potcar_path)                



class MagenticAnisotropySphereFlow:

    def __init__(self, workdir, npoints, kgrid, nbands, nodes, ppn, ref_orient, ldaul, magmom, Uparam, Jparam, encut, potcar_path, struct_path, ismear, name ="mae_calc", time_cl="12:00:00", time_ncl="01:40:00", sigma=0.2, nelect=None, cl_dir=None, angle_rng=None):

        """
        Computes the MCAE sphere for a defined structure.
        """

        self._workdir = workdir
        self._npoints = npoints
        self._name = name
        self._structure_path = struct_path
        self._potcar_path =  potcar_path
        self._reference_orientation = ref_orient
        self._collinear_calc = None
        self._cl_dir = cl_dir
        self._angle_rng = angle_rng or [[0, 2*np.pi],[0,np.pi]]         
 
        if self._cl_dir:
           pass
        else:
           self._cl_dir = self._workdir+"/"+"scf_cl"
        self._non_collinear_calcs = []

        def generate_spin_axes_h(npoints):
            golden_angle = (3 - np.sqrt(5)) * np.pi
            theta = golden_angle * np.arange(npoints)
            Sz = np.linspace(1/npoints-1, 1-1/npoints, npoints)
            radius = np.sqrt(1 - Sz * Sz)
            Sy = radius * np.sin(theta)
            Sx = radius * np.cos(theta)
            saxes_temp = np.round(np.array([Sx,Sy,Sz]).T, 4)
            saxes = saxes_temp.tolist()
            return saxes

        def generate_semi_spin_axes_h(npoints, thet_rng, phi_rng):
            num = math.floor(math.sqrt(npoints))
            thet = np.linspace(thet_rng[0], thet_rng[1], num)
            phi = np.linspace(phi_rng[0], phi_rng[1], num)
            u, v = np.meshgrid(thet, phi)
            x=np.cos(u)*np.sin(v)
            y=np.sin(u)*np.sin(v)
            z=np.cos(v)
            pts = np.array([x.flatten(), y.flatten(), z.flatten()])
            pts = pts.T 
            return pts.tolist() 

        if angle_rng: 
           self._saxes =  generate_semi_spin_axes_h(self._npoints, angle_rng[0], angle_rng[1])
        else: 
           self._saxes = generate_spin_axes_h(self._npoints)

	#self._saxes.append(self._reference_orientation)

        def set_calculations_h(cl_calc, cl_dir, saxes, workdir, nodes, ppn, time_cl, time_ncl, ismear, sigma, kgrid, encut, magmom, ldaul, Uparam, Jparam, nbands, nelect):
            if cl_calc:
                collinear_calc = None
            else:
               cl_settings = DefaultMagCLParameters(encut=encut, magmom=magmom, ldaul=ldaul, Uparam=Uparam, Jparam=Jparam)
               cl_settings.update_parallel_settings("flnm ", "run_cl.sh")
               cl_settings.update_parallel_settings("job_name", "cl_run")
               cl_settings.update_parallel_settings("nodes", nodes)
               cl_settings.update_parallel_settings("ppn", ppn)
               cl_settings.update_parallel_settings("max_time", time_cl)
               cl_settings.update_electronic_settings("ISMEAR", ismear)
               cl_settings.update_electronic_settings("SIGMA", sigma)
               cl_settings.update_electronic_settings("EDIFF", 1.0E-6)
               if nelect:
                    cl_settings.update_start_settings("NELECT", nelect)
               else:
                 pass

               collinear_calc = SCFCalculation(cl_dir, pseudo_par=None, kgrid=kgrid, name="scfcl", input_parameters=cl_settings)
            itr = 0
            non_collinear_calcs = []
            for spin_axis in saxes:
                ncl_settings = DefaultMagNCLParameters(encut=encut, spinaxis=spin_axis, ldaul=ldaul, Uparam=Uparam, Jparam=Jparam)
                ncl_settings.update_start_settings("NBANDS", nbands)
                ncl_settings.update_start_settings("LWAVE", ".FALSE.")
                ncl_settings.update_parallel_settings("flnm ", "run_ncl.sh")
                ncl_settings.update_parallel_settings("job_name", "ncl_run_"+str(itr))
                ncl_settings.update_parallel_settings("nodes", nodes)
                ncl_settings.update_parallel_settings("ppn", ppn)
                ncl_settings.update_parallel_settings("max_time", time_ncl)
                ncl_settings.update_parallel_settings("KPAR", None)
                ncl_settings.update_parallel_settings("exec", "vasp_ncl")
                ncl_settings.update_electronic_settings("ISMEAR", ismear)
                ncl_settings.update_electronic_settings("SIGMA", sigma)
                ncl_settings.update_electronic_settings("EDIFF", 1.0E-4)
                if nelect:
                     ncl_settings.update_start_settings("NELECT", nelect)
                else:
                  pass

                ncl_dir = workdir+"/"+"scf_ncl"+"/"+"scf_ncl_"+str(itr)
                ncl_calc = SCFCalculation(ncl_dir, pseudo_par=None, kgrid=kgrid, name="scfncl"+str(itr), input_parameters=ncl_settings)
                non_collinear_calcs.append(ncl_calc)
                itr += 1
            return [collinear_calc, non_collinear_calcs]


        [self._collinear_calc, self._non_collinear_calcs] = set_calculations_h(self._collinear_calc, self._cl_dir, self._saxes, self._workdir, nodes, ppn, time_cl, time_ncl, ismear, sigma, kgrid, encut, magmom, ldaul, Uparam, Jparam, nbands, nelect)

    def make_calculations(self):
        os.mkdir(self._workdir)
        self._collinear_calc.make_calculation(struct_path=self._structure_path, potcar_path=self._potcar_path)
        os.mkdir(self._workdir+"/"+"scf_ncl")
        for calc in self._non_collinear_calcs: calc.make_calculation(struct_path=self._structure_path, potcar_path=self._potcar_path)

    # def run_cl_calculation(self):
    #     self._collinear_calc.run_calculation()

    # def run_ncl_calculations(self):
    #     if self._collinear_calc._run_status != "finished":
    #         self._collinear_calc.update_run_status()
    #         self.run_ncl_calculations()
    #     else:
    #        cl_wvcr = self._collinear_calc._workdir+"/"+"WAVECAR"
    #        cl_chcr = self._collinear_calc._workdir+"/"+"CHGCAR"
    #        for calc in self._non_collinear_calcs:
    #            copyfile(cl_wvcr, calc._workdir)
    #            copyfile(cl_chcr, calc._workdir)
    #            calc.run_calculation()

    # def retrieve_mae_data(self):
    #     mae_data = []
    #     for calc in self._non_collinear_calcs:
    #         if calc._run_status != "finished":
    #            self._collinear_calc.update_run_status()
    #            self.retrieve_mae_data()
    #         else:
    #             spin = calc._input_settings._magnetic_settings["SAXIS"]
    #             energ = calc.get_total_energy()
    #             os.system("rm -r "+calc._workdir+"/"+"WAVECAR")
    #             os.system("rm -r "+calc._workdir+"/"+"CHGCAR")
    #             mae_data.append([spin[0], spin[1], spin[2], energ])
    #     write_maefile(mae_data)



	# class BandCalculation(Calculation):
	#         def __init__(self, charge_option, wfn,  istart, work_dir, structure, k_dim, pseudo_lists, k_path = None, nbnds = None, lorbit = True , smear = False , sigma = 0.01, isym = 0):
	#             """OUTPUT List:
	#                Istart:
	#                   - 0 = begin from scratch
	#                   - 1 = continue job w energy cutoff
	#                   - 2 = continue, restart w constant basis
	#                ICHARG:
	#                    -   0 (computes from initial wfn)
	#                    - = 1 extrapolate from old positions, reads from CHCAR
	#
	#               LORBIT =  10 ==> not decomposed DOS
	#               LORBIT = 11 ==> decomposed DOS
	#             """
	#             super(SCFCalculation, self).__init__(self, work_dir, structure, k_dim, pseudo_lists)
	#
	#             self._bands_options_ = [k_path, nbnds, charge_option]
	#
	#             self._scf_setting_ = [charge_option, prec, encut, nstep, epsilon, pseudo, n_elect.structure, smear, sigma, isym]
	#             self._ionic_  = []
	#             self._dos_ = [lorbit, dos_pts]
	#             self._bands_options_ = [k_path, nbnds, charge_option]
	#             self._magentic_options_ = [spin_polar, magmom, spin_orbit]
	#             self._hubbard_ = [dft_u, dudarev, ldaul, u_param, j_param, lda_mix]
	#             self._hybrid_ = [hf_calc, hf_fft, hf_max, hf_scren, hf_xch]
	#             self._rho_decomp_ = [par_chg, en_rng, ref_ef, over_k, over_bnd]
	#
	#         def make_calculation(self):
	#                 """Sets input parameters for ground state calculation and makes files for calculation"""
	#             temp_dir = os.getcwd()
	#             os.mkdir(self._workdir)
	#             print("Work Directory now in:" + self._workdir)
	#             os.chdir(self._workdir)
	#             make_incar(self._pseudos, self._input)
	#             make_poscar(self. )
	#             make_potcar(self._pseudos)
	#             make_runscript(self._parallelization)
	#             os.chdir(temp_dir)
	#


	#class PhononCalculation(Calculation):
	#    def __init__(self, input):
	#        """""
	#        Sets input for relaxation calculation
	#        """
	#        self._algorithm_ =
	#        self._bands_ =
	#        self._rho_   =
	#        self._spin_orbit_ =
	#        self._charge_compensate_ =
	#        self._spin_polar_  =
	#        super(RelaxationCalculation, self).__init__(self, work_dir, structure, k_dim, pseudo_list )



	#class EFieldCalculation(Calculation):


	# class SerialComputeFlow():
	# """
	# Sets up a serial set of compuations which iterate over a desired degree of freedom
	# (i.e. strain, dopant type/position, substrate)
	# """
	#     def __init__(self, workdir, dir_prefix="scf_calc",name="series_calculation" ):
	#     """
	#     Sets up a serial set of compuations which iterate over a desired degree of freedom
	#     (i.e. strain, dopant type/position, substrate)
	#     """
	#         self._workdir = workdir
	#         self._dir_prefix = dir_prefix
	#         self._dir_names = []
	#         self._ncalc = 0
	#         self._param_list = []
	#         self._calc_list = []

	    # def setup_calc_series(self):
	    #     os.chdir(self._workdir)
	    #     param_list = self._param_list
	    #     for i in range(self._ncalc):
	    #         dir_name = self._dir_prefix + param_list[i]
	    #         self._dir_names.append(dirname)
	    #         os.mkdir(dirname)

	    # def make_series(self):
	    #     if len(self._calc_list) == 0:
	    #         pass
	    #     else:
	    #         for calc in self._calc_list:
	    #             calc.make_calculation()
	    # def add_calc(self, param):
	    #
	    # def remove_calc(self):
	    #
	    # def run_status(self):
	    #
	    #
	    #

	# class ConvergeTest(SerialComputeFlow):
	#     """
	#     Sets up a serial set of compuations which iterate over a desired degree of freedom
	#     (i.e. strain, dopant type/position, substrate)
	#     """

	    # """OUTPUT List:
	    #    - CHG: chrge density, lattice vecctors, coords.
	    #    - DOSCAR: DOS
	    #    - EIGENVAL: bands
	    #    - IBZKPT: BZ
	    #    - LOCPOT: local potential
	    #    - OSZICAR: information at each nstep
	    #    - OUTCAR: outputfile, main
	    #    - PARCHG: partial charage density
	    #    - PROCAR: site projected wfn cahracter
	    #    - WAVECAR: wavefunctions and coefficients, eigenvalues
	    # """
	    #
	    #     self._algorithm_ =
	    #     self._bands_ =
	    #     self._rho_   =
	    #     self._spin_orbit_ =
	    #     self._charge_compensate_ =
	    #     self._spin_polar_  =
	    #     self._dos_ =
	    #     super(SerialComputeFlow, self).__init__(inialized inputs, )

	#class TimingTest(SerialComputeFlow):#
	#    """
	    # Sets up a serial set of compuations which iterate over a desired degree of freedom
	    # (i.e. strain, dopant type/position, substrate)
	    # """
	    #     self._algorithm_ =
	    #     self._bands_ =
	    #     self._rho_   =
	    #     self._spin_orbit_ =
	    #     self._charge_compensate_ =
	    #     self._spin_polar_  =
	    #     self._dos_ =
	    #     super(SerialComputeFlow, self).__init__(inialized inputs, )
