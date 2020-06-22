import numpy as np
import pymatgen as pg
import scipy as sp
import os
from io import *

class Calculation():
    def __init__(self, work_dir, structure, k_dim, pseudo_list):
        """
        Calculation superclass defines and produces the sequence of input files and directories
        needed to compute a specified property in VASP. Input parameters, structure, pseudopotentials,
        and degrees of freedom initially specified.
        """
        try:
          os.mkdir(work_dir)
        except OSError as error:
          print("Directory not found, try another path.")

        self._workdir_ = work_dir
        self._structure_ = structure
        self._kgrid_ = k_dim
        self._pseduo_dir_ = '/global/common/sw/cray/cnl7/haswell/vasp/pseudopotentials/PBE/potpaw_PBE'
        self._pseudos_ = pseudolist
        self._parallelization_ = ['knl','regular', 0,0, k_dim[0]*k_dim[1]*k_dim[2]] ##[computer, calc_type, nodes, nthreads, ppn]

    def set_structure(self, structure):
        """Initializes structure for calculation"""
        self._structure_ = structure

    def set_pseudos(self, pseudo_dir, pseudo_list): ##: specify criteria definining pseudo potential
        """Sets input parameters for ground state calculation"""
        self._pseduo_dir_ = pseudo_dir
        self._pseudos_ = pseudolist

    def set_parallelization(self,  calc_type, machine, nodes, nthreads, ppn): # **: specify criteria defining paralleliation
        """Sets input parameters for ground state calculation"""
        assert calc_type == 'regular' or calc_type == 'colinear' or calc_type == 'noncolinear', 'Calcualtion type must be regular, colinear, or noncolinear'
        self._parallelization_ = ['knl', calc_type, nodes, nthreads, ppn]

class RelaxationCalculation(Calculation):
    def __init__(self, relax_type, encut, nlec, nstep, work_dir, structure, k_dim, pseudo_list, epsilon=0.0001, nsteps=60, smear=0):
        """
        RelaxationCalculation subclass defines the state of a relaxation calculation executed in VASP, provided the structure, potentials, and kmesh
        defined in the superclass.  RelaxationCalculation is a subclass of Calculation. Calculation type can be an ionic relaxation, volume relaxation,
        or combination of the two. Optimization excecuted via conjugate-gradient algorithm.
        """
        assert relax_type == 'v' or relax_type == 'p' or relax_type == 's' or relax_type == 'vp' or relax_type == 'vs' or relax_type == 'ps' or relax_type == 'vps', 'relax_type must positions, volume, shape, or any combination of the 3'
        self._relax_type_ = relax_type
        self._encut_ = encut
        self._nsteps_ = nsteps
        self._epsilon_ = epsilon
        self._smear_ = smear
        super(RelaxationCalculation, self).__init__(self, work_dir, structure, k_dim, pseudo_list)
        self._num_elec_ = structure._num_elec_


    def set_calculation(self):
        """Sets input parameters for ground state calculation and makes files for calculation"""
        temp_dir = os.getcwd()
        os.mkdir(self._workdir_)
        print("Work Directory now in:" + self._workdir_)
        os.chdir(self._workdir_)
        make_incar(self._structure_, self._pseudos_, self._input_)
        make_poscar(self._structure_)
        make_potcar(self._pseudos_)
        make_runscript(self._parallelization_)
        os.chdir(temp_dir)

#    def get_relaxedstructure(self, out_dir = self._workdir_):
#        """Retrieves relaxed structure in the form of a POSCAR file"""
#        os.chdir(self._workdir_)
#        output = open("OUTCAR", "r")


class SCFCalculation(Calculation):
    def __init__(self, input, bands = True, rho = True, spin_orbit = True, spin_polar = 'regular', spin_axis = [0,0,0]):
        """""
        Sets input for relaxation calculation
        """
        self._bands_ = bands
        self._rho_   = rho
        self._spin_orbit_ = True
        self._spin_polar_  = 'regular'
        self._dos_ =
        super(RelaxationCalculation, self).__init__(bands,  )

        def set_calculation(self, ):


class PhononCalculation(Calculation):
    def __init__(self, input):
        """""
        Sets input for relaxation calculation
        """
        self._algorithm_ =
        self._bands_ =
        self._rho_   =
        self._spin_orbit_ =
        self._charge_compensate_ =
        self._spin_polar_  =
        self._dos_ =
        super(RelaxationCalculation, self).__init__(inialized inputs, )

class SerialComputeFlow():
"""
Sets up a serial set of compuations which iterate over a desired degree of freedom
(i.e. strain, dopant type/position, substrate)
"""
    def __init__(self, input):
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
        self._algorithm_ =
        self._bands_ =
        self._rho_   =
        self._spin_orbit_ =
        self._charge_compensate_ =
        self._spin_polar_  =
        self._dos_ =


class ConvergeTest(SerialComputeFlow):
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
        self._algorithm_ =
        self._bands_ =
        self._rho_   =
        self._spin_orbit_ =
        self._charge_compensate_ =
        self._spin_polar_  =
        self._dos_ =
        super(SerialComputeFlow, self).__init__(inialized inputs, )

class TimingTest(SerialComputeFlow):
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
        self._algorithm_ =
        self._bands_ =
        self._rho_   =
        self._spin_orbit_ =
        self._charge_compensate_ =
        self._spin_polar_  =
        self._dos_ =
        super(SerialComputeFlow, self).__init__(inialized inputs, )

class MagenticAnisotropyFlow(SerialComputeFlow):
    """
    Sets up a serial set of compuations which iterate over a desired degree of freedom
    (i.e. strain, dopant type/position, substrate)
    """
        self._algorithm_ =
        self._bands_ =
        self._rho_   =
        self._spin_orbit_ =
        self._charge_compensate_ =
        self._spin_polar_  =
        self._dos_ =
        super(RelaxationWorkflow, self).__init__(inialized inputs, )
