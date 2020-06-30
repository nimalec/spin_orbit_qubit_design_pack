import numpy as np
import pymatgen as pg
import scipy as sp

class Atom():

   def __init__(self):
       """ Class Atom defined by an atom's species, cartesian coordinate, number of electrons, and mass
       """
       self._species =
       self._coordinate =
       self._num_elec =
       self._mass =

   def translate_atom(self, trans_vector):
        """specify lattice vector of transation"""

   def cart2frac(self):
           """converts from cartesian to fractional coordinates"""

   def frac2cart(self):
           """converts from fractional to cartesian coordinates"""

class Basis():
   def __init__(self):
          """Class Basis defined by a set of Atoms in an irreducible representation"""
          self._atomlist =
          self._natoms =


class Lattice():
    def __init__(self):
        """Class Lattice defined by a set of """
          self._avect =
          self._bvect =
          self._cvect =
    def scale(self):

    def translate(self):

    def rotate(self):



class Supercell():
    def __init__(self):
        """Class Lattice defined by a set of Atoms in an irreducible representation"""


class Defect():
    def __init__(self):
        """Class Lattice defined by a set of Atoms in an irreducible representation"""


class Structure():
    def __init__(self):
        """Class Lattice defined by a set of Atoms in an irreducible representation"""
        self._species =
        self._formula =
        self._atoms =
        self._basis =
        self._cell_dimension =
        self._lattice  =
        self._basis =
        self._spacegroup =
        self._nelect =
        self._chargestate
        self._defects =
        self._polar_axes =
