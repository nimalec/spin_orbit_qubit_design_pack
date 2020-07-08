import numpy as np
import pymatgen as pg
import scipy as sp
import os
from io import *

class Calculation():
    def __init__(self, work_dir, structure, k_dim, pseudo_dir = '/global/common/sw/cray/cnl7/haswell/vasp/pseudopotentials/PBE/potpaw_PBE', pseudo_list = [     ], start_setting = [], parallelization = [group, partition, num_nodes, time_out]):
        """
        Calculation superclass defines and produces the sequence of input files and directories
        needed to compute a specified property in VASP. Input parameters, structure, pseudopotentials,
        and degrees of freedom initially specified.

        parallelization = [group, partition, num_nodes, ppn, time_out]
        """
        try:
          os.mkdir(work_dir)
        except OSError as error:
          print("Directory not found, try another path.")

        self._workdir = work_dir
        self._start_calc = start_setting
        self._structure = structure
        self._kgrid = k_dim
        self._pseudo_dir =
        self._pseudos = pseudolist
        self._parallelization = [jobname, machine, partition, nodes, ntasks, maxtime, ncore, kpar]

    def set_structure(self, structure):
        """Initializes structure for calculation"""
        self._structure_ = structure

    def set_pseudos(self, pseudo_dir, pseudo_list): ##: specify criteria definining pseudo potential
        """Sets input parameters for ground state calculation"""
        self._pseduo_dir = pseudo_dir
        self._pseudos = pseudolist

    def set_parallelization(self,  calc_type, machine, nodes, nthreads, ppn): # **: specify criteria defining paralleliation
        """Sets input parameters for ground state calculation"""
        assert calc_type == 'regular' or calc_type == 'colinear' or calc_type == 'noncolinear', 'Calcualtion type must be regular, colinear, or noncolinear'
        self._parallelization = ['knl', calc_type, nodes, nthreads, ppn]

    def get_run_time(self):


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
        self._num_elec = structure._num_elec_
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
    def __init__(self, charge_option, wfn,  istart, work_dir, structure, k_dim, pseudo_lists, k_path = None, nbnds = None, lorbit = True , smear = False , sigma = 0.01, isym = 0):
        """OUTPUT List:
           Istart:
              - 0 = begin from scratch
              - 1 = continue job w energy cutoff
              - 2 = continue, restart w constant basis
           ICHARG:
               -   0 (computes from initial wfn)
               - = 1 extrapolate from old positions, reads from CHCAR

          LORBIT =  10 ==> not decomposed DOS
          LORBIT = 11 ==> decomposed DOS
        """
        super(SCFCalculation, self).__init__(self, work_dir, structure, k_dim, pseudo_lists, )

        self._electronic_settings = [charge_option, prec, encut, nstep, epsilon, pseudo, n_elect.structure, smear, sigma, isym]
        self._ionic_settings  = []
        self._dos_settings = [lorbit, dos_pts]
        self._magentic_options_ = [spin_polar, magmom, spin_orbit]
        self._hubbard_ = [dft_u, dudarev, ldaul, u_param, j_param, lda_mix]
        self._hybrid_ = [hf_calc, hf_fft, hf_max, hf_scren, hf_xch]
        self._rho_decomp_ = [par_chg, en_rng, ref_ef, over_k, over_bnd]
        self._input_ = [self._electronic_setting_, self._parallelization_, self._ionic_ , self._dos_, self._magentic_options_, self._hubbard_, self._hybrid_, self._rho_decomp_]

        def make_calculation(self):
            """Sets input parameters for ground state calculation and makes files for calculation"""
            temp_dir = os.getcwd()
            os.mkdir(self._workdir)
            print("Work Directory now in:" + self._workdir_)
            os.chdir(self._workdir)
            make_incar(self._input_)
            make_poscar(self._structure_)
            make_potcar(self._pseudos_)
            make_runscript(self._parallelization_)
            os.chdir(temp_dir)

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
    def __init__(self, workdir, paramter,   ):
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
        self._workdir = workdir
        self._parameter = paramter
        self._dir_prefix = paramter
        self._dir_names = []
        self._ncalc = 0
        self._param_list =
        self._calc_list = []
        self._runscript =

    def setup_calc_series(self):
        os.chdir(self._workdir)
        param_list = self._param_list
        for i in range(self._ncalc):
            dir_name = self._dir_prefix + param_list[i]
            self._dir_names.append(dirname)
            os.mkdir(dirname)

    def add_calc(self):

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
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
    def __init__(self, cl_calc, polar_range, azimuth_range, ref_orient):
        self._reference_orientation = ref_orient
        self._cl_calc  =  cl_calc
        self._collinear_calc = cl_calc
        self._azimuth_range = polar_range
        self._polar_range = azimuth_range,
        self._symmetryreduce = [False, None, None]
        def generate_spin_axes_h(self):
            spin_list = []
            itr = 0
            for i in self._azimuth_range:
                for j in  self._polar_range:
                    spin_list.append([itr, np.sin(i)*np.cos(j), np.sin(i)*np.sin(j), cos(j)])
                    itr += 1
            return spin_list
        self._spin_list = self.generate_spin_axes_h()
        super(SerialComputeFlow, self).__init__(inialized inputs, reference_orient)
           
