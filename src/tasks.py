import numpy as np
#import pymatgen as pg
#import scipy as sp
import os
#from io import *

__author__ = 'Nima Leclerc'
__email__ = 'nleclerc@lbl.gov'

# class Calculation():
#     def __init__(self, calc_name, work_dir, structure, k_dim, pseudo_par , parallel = ["run_cl.sh", "nano", "etna", 24, 24, 24, 8, 4]):
#         """
#         Calculation superclass defines and produces the sequence of input files and directories
#         needed to compute a specified property in VASP. Input parameters, structure, pseudopotentials,
#         and degrees of freedom initially specified.
#         ** parallel = [jobname, machine, partition, nodes, ppn, maxtime, ncore, kpar]
#         **
#         **
#
#         """
#         assert type(calc_name) is str, "calc_name not a String!"
#         assert type(work_dir) is str, "work_dir not a String!"
#         assert type(structure) is Structure, "structure not a Structure!"
#         assert len(k_dim) == 3, "k_dim must have just 3 elements!"
#         assert len(pseudo_par) == 2, "pseudo_par must have 2 elements!"
#         assert len(parallel) == 8, "parallel must have 8 elements!"
#
#         try:
#           os.mkdir(work_dir)
#         except OSError as error:
#           print("Directory not found, try another path.")
#
#         self._calc_settings = {"name": calc_name, "directory": work_dir}
#         self._structure = struct
#         self._pseudos = pseudo_par
#         self._parallelization = parallel
#
#     def set_structure(self, structure):
#         """Initializes structure for calculation"""
#         self._structure_ = structure

class InputParameters():
        def __init__(self, name="input_param", start_settings=None,
        parallel_settings=None, electronic_settings=None, magnetic_settings=None,
        hybrid_settings=None, hubbard_settings=None, misc_settings=None):
            """
            Sets standard input parameters for a VASP calculation.

            **Args:

            name (str): name of input paramter settings [default = "input_param"]
            start_settings (dict): start settings for calculation [default = ]
            parallel_settings (dict): parallization settings for calculation [default = ]
            electronic_settings (dict): electronic structure settings for calculation [default = ]
            ionic_settings (dict): ionic settings for calculation, used for relaxations and structure optimization [default=None]
            magnetic_settings (dict): magnetic structure settings for calculation [default=None]
            hybrid_settings (dict): hybrid/hse settings for accurate band calculation [default=None]
            hubbard_settings (dict): hubbard calculation settings for localized d/f orbital predicitons [default=None]
            misc_settings_settings (dict): miscalaneous settings for calculation, can be any VASP setting [default=None]

            """

            self._name = name
            self._start_settings = start_settings or {"nwrite": , "istart": , "iniwav": , "icharg": , "nelect": , "lorbit": , "nedos": , "loptics": , "lelf":, "lvhar": , "rwigs": , "lvtof": }
            self._parallel_settings = parallel_settings or {"flnm": , "job_name": , "machine":, "partition":, "nodes":  ,"ppn": , "ppn": , "max_time": , "ncore": , "kpar": }
            self._electronic_settings = electronic_settings or  {"rec_level": , "algo":, "encut": ,"nelm": ,"nelmin":, "ediff":, "sigma": ,"lasph": , "lreal":, "addgrid": , "bmaxmix": , "bmix": }
            self._ionic_settings = ionic_settings
            self._magnetic_settings = magnetic_settings
            self._hybrid_settings = hybrid_settings
            self._hubbard_settings = hubbard_settings
            self._misc_settings = misc_settings

        def get_input_settings(self):
            input_settings = {"name": name, "start": self._start_settings, "parallel": self._parallel_settings ,
            "electronic": self._electronic_settings, "magnetic": self._magnetic_settings, "hybrid": self._hybrid_settings, "hubbard": self._hubbard_settings}
            return input_settings

        def update_start_settings(self, key, value):
            if key in self._start_settings:
                self._start_settings[key] = value
                print "key" + "key changed to" + str(value)
            else:
                print "key does not exist!! keys include: {charge_option, prec, encut, nstep, epsilon, pseudo, n_elect.structure, smear, sigma, isym}"

        def update_parallel_settings(self, key, value):
            if key in self._parallel_settings:
                self._parallel_settings[key] = value
                print "key" + "key changed to" + str(value)
            else:
                print "key does not exist!! keys include: {flnm , job_name , machine, partition, nodes  ,ppn, max_time , ncore,  kpar}"

        def update_electronic_settings(self, key, value):
            if key in self._electronic_settings:
                self._electronic_settings[key] = value
                print "key" + "key changed to" + str(value)
            else:
                print "key does not exist!! keys include: {prec_level, algo, encut , nelm,nelmin, ediff, sigma, lasph, lreal, addgrid, bmaxmix, bmix}"

        def update_ionic_settings(self, key, value):
            if self._ionic_settings:
              if key in self._ionic_settings:
                self._ionic_settings[key] = value
                print "key" + "key changed to" + str(value)
              else:
                print "key does not exist!! keys include: {ediff ,nsw, ibrion ,isif, isym, nblock,  kblock}"
            else:
              print "magnetic settings not present!"

        def update_magnetic_settings(self, key, value):
            if self._magnetic_settings:
              if key in self._magnetic_settings:
                self._magnetic_settings[key] = value
                print "key" + "key changed to" + str(value)
              else:
                print "key does not exist!! keys include: {ispin, magmom, nupdown, saxis, lsorbit,noncollinear}"
            else:
              print "magnetic settings not present!"

        def update_hybrid_settings(self, key, value):
            if self._hybrid_settings:
              if key in self._hybrid_settings:
                  self._hybrid_settings[key] = value
                  print "key" + "key changed to" + str(value)
              else:
                  print "key does not exist!! keys include: {lhfcalc, precfock, nkred, algo, time, hflmax, hfscreen, aexx}"
            else:
              print "hybrid settings not present!"

        def update_hubbard_settings(self, key, value):
            if self._hubbard_settings:
              if key in self._hubbard_settings:
                  self._hubbard_settings[key] = value
                  print "key" + "key changed to" + str(value)
              else:
                  print "key does not exist!! keys include: {ldau, ldatype, ldaul, dlauu, ldauj, lmaxmix}"
           else:
             print "hybrid settings not present!"


# class DefaultOptimizationParameters(InputParameters):
#         def __init__(self, ):


# class DefaultSCFParameters(InputParameters):
#         def __init__(self, ):

# class DefaultSCFUParameters(InputParameters):
#         def __init__(self, ):

# class DefaultSCFHSEParameters(InputParameters):
#         def __init__(self, ):

# class DefaultMagCLParameters(InputParameters):
#         def __init__(self, ):

# class DefaultMagNCLParameters(InputParameters):
#         def __init__(self, ):

# class DefaultBandsParameters(InputParameters):
#         def __init__(self, ):





class RelaxationCalculation(Calculation):
    def __init__(self, relax_type, encut, nlec, nstep, work_dir, structure, k_dim, pseudo_list, epsilon=0.0001, nsteps=60, smear=0):
        """
        RelaxationCalculation subclass defines the state of a relaxation calculation executed in VASP, provided the structure, potentials, and kmesh
        defined in the superclass.  RelaxationCalculation is a subclass of Calculation. Calculation type can be an ionic relaxation, volume relaxation,
        or combination of the two. Optimization excecuted via conjugate-gradient algorithm.
        """
        assert relax_type == 'v' or relax_type == 'p' or relax_type == 's' or relax_type == 'vp' or relax_type == 'vs' or relax_type == 'ps' or relax_type == 'vps', 'relax_type must positions, volume, shape, or any combination of the 3'
        super(RelaxationCalculation, self).__init__(self, work_dir, structure, k_dim, pseudo_list)
        self._relax_type = relax_type
        self._encut = encut
        self._nsteps = nsteps
        self._epsilon = epsilon
        self._smear = smear
        self._num_elec = structure._num_elec
        self._input = [self._relax_type, self._encut, self._nsteps, self._epsilon, self._smear, self._num_elec]

    def make_calculation(self):
        """Sets input parameters for ground state calculation and makes files for calculation"""
        temp_dir = os.getcwd()
        os.mkdir(self._workdir)
        print("Work Directory now in:" + self._workdir)
        os.chdir(self._workdir)
        make_incar(self._pseudos, self._input)
        make_poscar(self._structure)
        make_potcar(self._pseudos)
        make_runscript(self._parallelization)
        os.chdir(temp_dir)

class SCFCalculation(Calculation):
     def __init__(self, electronic_setings, calc_name, work_dir, structure, k_dim, pseudo_param, start_settings = [2, 1, 1, 1, 11],  ionic_settings = [None, 0, -1 ,2, -1, 1, 0] ,magnetic_settings = [1, [0,0,0],[0,0,1], False ,False], hubbard_settings = [False, 2, [],[], 4], hybrid_settings = [lhfcalc = False, "F" , 2, "damped", 0.3, 4, .207, .25] , rho_settings = [True, [-10, 0], -3, 1, 20] , parallel = ["run_cl.sh", "nano", "etna", 24, 24, 24, 8, 4]):
        """OUTPUT List:
        start_settings = [nwrite ,istart ,iniwav ,icharg, nelect, lorbit, nedos, lvhar, rwigs, loptics, lvtof, lelf]
        electronic_settings = [prec_level, algo, encut , nelm,nelmin, ediff, sigma,   lasph, lreal, addgrid, bmaxmix, bmix]
        ionic_settings = [ediff, nsw, ibrion ,isif, isym, nblock, kblock]
        magnetic_settings = [ispin, magmom, nupdown, saxis, lsorbit,noncollinear]
        hubbard_settings = [ldau, ldatype, ldaul, dlauu, ldauj, lmaxmix]
        hybrid_settings = [lhfcalc, precfock, nkred, algo, time, hflmax, hfscreen, aexx]
        rho_decomp = [lpard, eint, dbmod, kpuse, iband]
        """
        super(SCFCalculation, self).__init__(self, calc_name, work_dir, structure, k_dim, pseudo_par, parallel)

        # assert type( ) is
        # assert type( ) is
        # assert type( ) is
        # assert type(  ) is
        # assert type( ) is

        self._start_settings = start_settings
        self._electronic_settings = electronic_settings
        self._ionic_settings  = ionic_settings
        self._magentic_settings = magnetic_setting
        self._hubbard_settings = hubbard_settings
        self._hybrid_settings = hybrid_settings
        self._rho_settings = rho_settings

        def make_calculation(self):
            """Sets input parameters for ground state calculation and makes files for calculation"""
            os.mkdir(self._workdir)
            print("Work Directory now in:" + self._workdir_)
            make_incar(self._workdir, self._input_settings)
            make_poscar(self._workdir, self._structure)
            make_potcar(self._workdir, self._pseudos)
            make_runscript(self._workdir, self._parallelization)

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

        #def get_mag_moment(self):
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


class SerialComputeFlow():
"""
Sets up a serial set of compuations which iterate over a desired degree of freedom
(i.e. strain, dopant type/position, substrate)
"""
    def __init__(self, workdir, paramter  ):
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
        self._workdir = workdir
        self._parameter = parameter
        self._dir_prefix = paramter
        self._dir_names = []
        self._ncalc = 0
        self._param_list = []
        self._calc_list = []
        self._runscript = make_serial_runscript()

    def setup_calc_series(self):
        os.chdir(self._workdir)
        param_list = self._param_list
        for i in range(self._ncalc):
            dir_name = self._dir_prefix + param_list[i]
            self._dir_names.append(dirname)
            os.mkdir(dirname)

    def make_series(self):
        if len(self._calc_list) == 0:
            pass
        else:
            for calc in self._calc_list:
                calc.make_calculation()



    def add_calc(self, para):

    def remove_calc(self):

    def run_status(self):




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

class MagenticAnisotropyFlow(SerialComputeFlow):

    def __init__(self, cl_calc, polar_range, azimuth_range, ref_orient):
        """
        Sets up a serial set of compuations which iterate over a desired degree of freedom
        (i.e. strain, dopant type/position, substrate)
        """
        self._reference_orientation = ref_orient
        self._collinear_calc = cl_calc
        self._azimuth_range = polar_range
        self._polar_range = azimuth_range,
        self._symmetryreduce = [False, None, None]

        def generate_spin_axes_h(self):
            spin_list = []
            itr = 0
            for phi in self._azimuth_range:
                for theta in  self._polar_range:
                    spin_list.append([itr, np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])
                    self._param_list.append([itr, np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(phi)])
                    itr += 1
            return spin_list
        self._spin_list = self.generate_spin_axes_h()

        def set_calcualtions_h(self):
            calc_list = []
            for spin in self._spin_list:
                temp_magnetic_settings = ['ncl', spin, None]
                temp_calc = SCFCalculation(self._collinear_calc.start_settings, self._collinear_calc.electronic_settings,self._collinear_calc._ionic_settings,temp_magnetic_settings, self._collinear_calc._hubbard_settings , self._collinear_calc._hybrid_settings,  self._collinear_calc._rho_settings)
                calc_list.append(temp_calc)
                self._ncalc += 1
            self._calc_list = calc_list
        super(SerialComputeFlow, self).__init__(inialized inputs, reference_orient)
        self.set_calcualtions_h()

    def add_spin(self, spin):
        """
        Add spin to spin list and make directory
        """
        temp_magnetic_settings = ['ncl', spin, None]
        temp_calc = SCFCalculation(self._collinear_calc.algorithm_settings, self._collinear_calc.electronic_settings, temp_magnetic_settings, self._collinear_calc.hubbard_settings)
        calc_list.append(temp_calc)
        self._ncalc += 1