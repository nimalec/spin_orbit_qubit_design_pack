import numpy as np
import pymatgen as pg
import scipy as sp
from io import makepotcar, makeposcar, makeincar

class Calculation():
    def __init__(self, work_dir, structure, k_dim, pseudo_list):
        """
        Workflow class defines and produces the sequence of input files and directories
        needed to compute a specified property. Input parameters, structure, pseudopotentials,
        and degrees of freedom initially specified.
        """
        self._workdir_ = work_dir
        self._structure_ = structure
        self._kgrid_ = k_dim
        self._pseduo_dir_ = '/global/common/sw/cray/cnl7/haswell/vasp/pseudopotentials/PBE/potpaw_PBE'
        self._pseudos_ = pseadolist
        self._parallelization_ = ['knl','regular', 0,0, k_dim[0]*k_dim[1]*k_dim[2]] ##[computer, calc_type, nodes, nthreads, ppn]

    def set_structure(self, structure):
        """Initializes structure for calculation"""
        self._structure_ = structure

    def set_pseudos(self, pseudo_dir, pseudo_list):
        """Sets input parameters for ground state calculation"""
        self._pseduo_dir_ = pseudo_dir
        self._pseudos_ = pseadolist

    def set_parallelization(self,  calc_type, machine, nodes, nthreads, ppn):
        """Sets input parameters for ground state calculation"""
        assert calc_type == 'regular' or calc_type == 'colinear' or calc_type == 'noncolinear', 'Calcualtion type must be regular, colinear, or noncolinear'
        self._parallelization_ = ['knl', calc_type, nodes, nthreads, ppn]

class RelaxationCalculation(Calculation):
    def __init__(self, relax_type, Ecut, nlec, nstep, work_dir, structure, k_dim, pseudo_list)
     """Sets input parameters for ground state calculation"""
        self._relax_type_ = 'regular'
        self._encut_ = Ecut
        self._epsilon_ = 0.0001
        self._num_elec_ = nelec
      super(RelaxationCalculation, self).__init__(nlec, nstep, epsilon, work_dir, structure, k_dim, pseudo_list)

    def set_calculation(self):
        """Sets input parameters for ground state calculation"""
        os.makedir()
        os.chddir(self._workdir_)
        make_incar(self._structure_, self._pseudos_, self._input_)
        make_poscar(self._structure_)
        make_potcar(self._pseudos_)

    def run_calculation(self):
        """Runs relaxation calcualtion in directory"""


    def get_relaxedstructure(self):
        """Retrieves relaxed structure in the form of a POSCAR"""



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


class ElectronPhononCalculation(Calculation):
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
        super(RelaxationWorkflow, self).__init__(inialized inputs, )

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
