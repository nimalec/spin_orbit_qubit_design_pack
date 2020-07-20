import numpy as np
#import pymatgen as pg
#import scipy as sp
import os
#from timer import Timer
#import time
from in_out import make_runscript_h, make_incar_h

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
            self._start_settings = start_settings or {"nwrite": 2, "istart": 1, "iniwav": 1,
             "icharg": None, "nelect": None, "lorbit": 11,
              "nedos": 1000, "loptics": False, "lelf": None, "lvhar": None, "rwigs": None, "lvtof": None}
            self._parallel_settings = parallel_settings or {"flnm": "run_scf.sh", "job_name": "scf_std", "machine": "nano" ,
             "partition": "etna", "nodes": 4,"ppn": 24,
              "max_time": "24:00:00", "ncore": 8, "kpar": 2, "exec": "vasp_std"}
            self._electronic_settings = electronic_settings or  {"algo": "Normal", "encut": 800,
            "nelm": 200, "nelmin": 4, "ediff": 10E-5, "ismear": 1,
            "sigma": 0.2,"lasph": True, "lreal": "Auto", "addgrid": True, "maxmix": 100, "bmix": 1.5}
            self._ionic_settings = ionic_settings
            self._magnetic_settings = magnetic_settings
            self._hybrid_settings = hybrid_settings
            self._hubbard_settings = hubbard_settings
            self._misc_settings = misc_settings

        @classmethod
        def get_input_settings(self):
            input_settings = {"name": name,
             "start": self._start_settings, "parallel": self._parallel_settings ,
            "electronic": self._electronic_settings, "magnetic": self._magnetic_settings,
            "hybrid": self._hybrid_settings, "hubbard": self._hubbard_settings, "misc_setting": self._misc_settings}
            return input_settings

        @classmethod
        def update_start_settings(self, key, value):
            if key in self._start_settings:
                self._start_settings[key] = value
                print("key" + "key changed to" + str(value))
            else:
                print("key does not exist!! keys include: {charge_option, prec, encut, nstep, epsilon, pseudo, n_elect.structure, smear, sigma, isym}")

        @classmethod
        def update_parallel_settings(self, key, value):
            if key in self._parallel_settings:
                self._parallel_settings[key] = value
                print("key" + "key changed to" + str(value))
            else:
                print("key does not exist!! keys include: {flnm , job_name , machine, partition, nodes  ,ppn, max_time , ncore,  kpar}")

        @classmethod
        def update_electronic_settings(self, key, value):
            if key in self._electronic_settings:
                self._electronic_settings[key] = value
                print("key" + "key changed to" + str(value))
            else:
                print("key does not exist!! keys include: {prec_level, algo, encut , nelm,nelmin, ediff, sigma, lasph, lreal, addgrid, bmaxmix, bmix}")

        @classmethod
        def update_ionic_settings(self, key, value):
            if self._ionic_settings:
              if key in self._ionic_settings:
                self._ionic_settings[key] = value
                print("key" + "key changed to" + str(value))
              else:
                print("key does not exist!! keys include: {ediff ,nsw, ibrion ,isif, isym, nblock,  kblock}")
            else:
              print("magnetic settings not present!")

        @classmethod
        def update_magnetic_settings(self, key, value):
            if self._magnetic_settings:
              if key in self._magnetic_settings:
                self._magnetic_settings[key] = value
                print("key" + "key changed to" + str(value))
              else:
                print("key does not exist!! keys include: {ispin, magmom, nupdown, saxis, lsorbit,noncollinear}")
            else:
              print("magnetic settings not present!")

        @classmethod
        def update_hybrid_settings(self, key, value):
            if self._hybrid_settings:
              if key in self._hybrid_settings:
                  self._hybrid_settings[key] = value
                  print("key" + "key changed to" + str(value))
              else:
                  print("key does not exist!! keys include: {lhfcalc, precfock, nkred, algo, time, hflmax, hfscreen, aexx}")
            else:
              print("hybrid settings not present!")

        @classmethod
        def update_hubbard_settings(self, key, value):
            if self._hubbard_settings:
              if key in self._hubbard_settings:
                  self._hubbard_settings[key] = value
                  print("key" + "key changed to" + str(value))
              else:
                  print("key does not exist!! keys include: {ldau, ldatype, ldaul, dlauu, ldauj, lmaxmix}")
            else:
             print("hybrid settings not present!")


class DefaultOptimizationParameters(InputParameters):
        def __init__(self, encut, name='relax_settings'):
            """
            Sets standard input parameters for a VASP calculation.

            **Args:

            encut (flt): planewave energy cutoff for calculation
            name (str): name for relaxation setting

            """

            ionic = {"ediff": 1E-17, "nsw": 20, "ibrion": 2,"isif": 2, "isym": -1, "nblock": 1,  "kblock": 20}
            InputParameters.__init__(self, ionic_settings=ionic, name=name)
            self.update_electronic_settings("encut", encut)

class DefaultSCFParameters(InputParameters):
         def __init__(self, encut, name='scf_settings'):
             """
             Sets standard input parameters for a VASP calculation.

             **Args:
               encut (flt): planewave energy cutoff for calculation
               name (str): name for relaxation setting

             """
             InputParameters.__init__(self, name=name)
             self.update_electronic_settings("encut", encut)

class DefaultSCFUParameters(InputParameters):
         def __init__(self, encut, ldaul, Uparam, Jparam, name='DFTU_settings'):
             """
             Sets standard input parameters for a VASP calculation.

             encut (flt): planewave energy cutoff for calculation
             name (str): name for relaxation setting

             """

             dftu_settings = {"ldau": Uparam, "ldatype": 2, "ldaul": ldaul, "ldauj": Jparam , "lmaxmix": 4}
             InputParameters.__init__(self, name=name, hubbard_settings=dftu_settings)
             self.update_electronic_settings("encut", encut)

# class DefaultSCFHSEParameters(InputParameters):
#         def __init__(self):
#             pass


class DefaultMagCLParameters(InputParameters):
         def __init__(self, encut, magmom, ldaul, Uparam, Jparam, name='DFTCL_settings'):
             """
             Sets standard input parameters for a VASP calculation.

             encut (flt): planewave energy cutoff for calculation
             name (str): name for relaxation setting

             """

             cl_settings =  {"ispin": 2, "magmom": magmom, "saxis": None, "lsorbit": None, "noncollinear": None}
             dftu_settings = {"ldau": Uparam, "ldatype": 2, "ldaul": ldaul, "ldauj": Jparam , "lmaxmix": 4}
             InputParameters.__init__(self, name=name, magnetic_settings=cl_settings, hubbard_settings=dftu_settings)
             self.update_electronic_settings("encut", encut)

class DefaultMagNCLParameters(InputParameters):
         def __init__(self, encut, spinaxis, ldaul, Uparam, Jparam, name='DFTCL_settings'):
             """
             Sets standard input parameters for a VASP calculation.

             encut (flt): planewave energy cutoff for calculation
             name (str): name for relaxation setting

             """
             ncl_settings =  {"ispin": 2, "magmom": None, "saxis": True, "lsorbit": True, "noncollinear": True}
             dftu_settings = {"ldau": Uparam, "ldatype": 2, "ldaul": ldaul, "ldauj": Jparam , "lmaxmix": 4}
             InputParameters.__init__(self, name=name, magnetic_settings=ncl_settings, hubbard_settings=dftu_settings)
             self.update_electronic_settings("encut", encut)

class SCFCalculation():
     def __init__(self, structure, workdir, kgrid, pseudo_par, name="scf_calc", encut=600, input_parameters=None):
         """
         Sets standard input parameters for a VASP calculation.

         **Args:

         structure (Structure):
         workdir (str):
         kgrid (list):
         pseudo_par (dict):
         name (str):
         input_paramters (InputParameters):
         encut (flt):
         """
         self._name  = name
         self._workdir = workdir
         self._structure = structure
         self._kmesh = kgrid
         self._pseudo_par = pseudo_par
         self._input_settings = input_parameters or DefaultSCFParameters(encut=encut)
         self._run_status = "unstarted" #can be {unstarted, running, crashed, or complete}
         self._jobid = None
         self._cputime = None
         self._tot_energy = None
         self._fermi = None


# input_settings = DefaultSCFParameters(encut=300)
# work_dir = "/Users/nimalec/Documents/Confidential Work /griffin_summer_2020 [Confidential]/spin_orbit_qubit_design_pack"
# make_incar_h(work_dir, input_settings)
#


     def make_calculation(self, pseudo_path=None, struct_path=None, run_script_path=None):
         """Sets input parameters for ground state calculation and makes files for calculation"""
         os.mkdir(self._workdir)
         print("Work Directory now in: " + self._workdir)
         make_incar(self._workdir, self._input_settings)
         make_poscar(self._workdir, self._structure)
         make_potcar(self._workdir, self._pseudos)
         make_runscript(self._workdir, self._parallelization)


     # def run_calculation(self):
     #     os.system("sbatch")
     #     self._run_status = "started"
     #     self._jobid = os.system() ##retrieves job id

    # def update_calc_status(self):
    #  if self._jobid: ##if job has started
    #      def is_finished_h(self):
    #          ## looks at output file to check if finished
    #          if "General timing" is present in OUTCAR:
    #            return true
    #          else:
    #              return false
    #        if os.system(self._jobid) is present:
    #            self.self._run_status = "running"
    #        else:
    #          if is_finished_h() == true:
    #              self.self._run_status = "Finished"
    #          else:
    #              self._run_status= "Unfinished"

     def get_total_energy(self):
         energ_list = []
         fl_nm = self._workdir + 'OUTCAR'
         isfile = os.path.isfile(fl_nm)
         if isfile == False:
           print("OUTCAR file not present! try to re-run the calculation.")
            pass
         else:
           with open(fl_nm, 'r') as f:
             for line in f.readlines():
               if 'TOTEN' in line:
                 energ_list.append(line)
         tot_energ = energ_list[len(energ_list)-1]
         return float(tot_energ[30:40])

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


     def get_fermi(self):
         fl_nm = self._workdir + 'OUTCAR'
         isfile = os.path.isfile(fl_nm)
         if isfile == False:
             print("OUTCAR file not present! try to re-run the calculation.")
             pass
         else:
            with open(fl_nm, 'r') as f:
              for line in f.readlines():
                if 'E-fermi :' in line:
                    fermi_str = line
        return float(fermi_str[12:18])

    #  def start_calculation(self):
    #      def is_finished_h(self):
    #          ## looks at output file to check if finished
    #          if "General timing" is present in OUTCAR:
    #            return true
    #          else:
    #              return false
    #
    #      self.make_calculation()
    #      self.run_calculation()
    #      start = time.time()
    #      while self.s_finished_h() is True:
    #          calc_time = time.time() - start
    #          if calc_time < 30:
    #              pass
    #          else:
    #              if int(calc_time)%30 !=0:
    #                  pass
    #              else:
    #                  self.update_calc_status()
    #      self._cputime = self.get_run_time()
    #      self._tot_energy = self.get_total_energy()
    #      self._fermi = self.get_fermi()

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
#
#     def make_serial_runscript(self):


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

# class MagenticAnisotropyFlow(SerialComputeFlow):
#
#     def __init__(self, cl_calc, polar_range, azimuth_range, ref_orient):
#         """
#         Sets up a serial set of compuations which iterate over a desired degree of freedom
#         (i.e. strain, dopant type/position, substrate)
#         """
#         self._reference_orientation = ref_orient
#         self._collinear_calc = cl_calc
#         self._azimuth_range = polar_range
#         self._polar_range = azimuth_range,
#         self._symmetryreduce = [False, None, None]
#
#         def generate_spin_axes_h(self):
#             spin_list = []
#             itr = 0
#             for phi in self._azimuth_range:
#                 for theta in  self._polar_range:
#                     spin_list.append([itr, np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])
#                     self._param_list.append([itr, np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])
#                     itr += 1
#             return spin_list
#         self._spin_list = self.generate_spin_axes_h()
#
#         def set_calcualtions_h(self):
#             calc_list = []
#             for spin in self._spin_list:
#                 temp_magnetic_settings = ['ncl', spin, None]
#                 temp_calc = SCFCalculation(self._collinear_calc.start_settings, self._collinear_calc.electronic_settings,self._collinear_calc._ionic_settings,temp_magnetic_settings, self._collinear_calc._hubbard_settings , self._collinear_calc._hybrid_settings,  self._collinear_calc._rho_settings)
#                 calc_list.append(temp_calc)
#                 self._ncalc += 1
#             self._calc_list = calc_list
#         super(SerialComputeFlow, self).__init__(inialized inputs, reference_orient)
#         self.set_calcualtions_h()
#
#     def add_spin(self, spin):
#         """
#         Add spin to spin list and make directory
#         """
#         temp_magnetic_settings = ['ncl', spin, None]
#         temp_calc = SCFCalculation(self._collinear_calc.algorithm_settings, self._collinear_calc.electronic_settings, temp_magnetic_settings, self._collinear_calc.hubbard_settings)
#         calc_list.append(temp_calc)
#         self._ncalc += 1
