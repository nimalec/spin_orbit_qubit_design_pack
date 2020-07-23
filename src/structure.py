import numpy as np
import scipy as sp
import pymatgen as pg

__author__ = 'Nima Leclerc'
__email__ = 'nleclerc@lbl.gov'


class Site:
    def __init__(self, species, coord, site_idx = None, charge=None, mag_moment=np.array([0,0,0]), isdefect=False):
        """
        Initializes site objects specifying chemistry and position.

        **Args:

        species (str): element symbol, must be a valid element symbol in the periodic table
        coord (ndarray): fractional coordinate of atom in lattice site
        site_idx (int): site index
        charge (float or int): charge state of element  [default = maximum oxidation number of element]
        mag_moment (ndarray): orientation of magnetic moment [default = None]
        isdefect(bool): True if defect site in lattice [default = False]
        """
        
        assert isinstance(species, str), "species must be a string!!"

        try:
            elem = pg.core.periodic_table.Element(species)
        except ValueError:
            print(species+" "+"must be an element!!")

        assert isinstance(coord, np.ndarray), "coord must be an numpy array!"
        assert len(coord)==3, "coord must be length 3!!"

        if charge:
            isinstance(charge, int) or isinstance(charge, float), "charge must be a float or int!"
        else:
            pass

        assert isinstance(isdefect, bool)

        assert isinstance(mag_moment, np.ndarray), "magmom must be an numpy array!"
        assert len(mag_moment)==3, "magbetic moment be length 3!!"

        self._mag_moment = mag_moment
        self._site_idx = site_idx or 1
        self._species = species
        self._coord = coord
        elem = pg.core.periodic_table.Element(species)
        self._charge = charge or elem.max_oxidation_state
        self._isdefect = False

    @classmethod
    def translate(self, displacement):
        """
        Translates site with provided displacement.
        **Args:

        displacement (ndarray): translation vector
        """
        assert isinstance(displacement, np.ndarray), "must be an ndarray!"
        assert len(displacement)==3, "displacement must be length 3!"
        temp = self._coord
        self._coord = temp + displacement

    @classmethod
    def set_charge(self, charge):
        """
        Sets charge state of site.

        **Args:

        charge (int or float): charge state of site
        """
        assert isinstance(charge, int) or isinstance(charge, float), "charge must be a float or int"
        self._charge = charge

    @classmethod
    def set_mag_moment(self, moment):
        """
        Sets magnetic moment of site.

        **Args:

        moment (ndarray): magnetic moment of site
        """

        assert isinstance(moment, np.ndarray), "must be an ndarray!"
        assert len(moment)==3, "displacement must be length 3!"
        self._mag_moment = moment

    @classmethod
    def set_defect(self, defect):
        """
        Sets defect state of site.

        **Args:

        defect (bool): defect state of site, True if site is a defect
        """
        self._isdefect = defect

    @property
    def species(self):
        """Getter function for species of site"""
        return self._species

    @property
    def coord(self):
        """Getter function for fractional coordinates of site"""
        return self._coord

    @property
    def mag_moment(self):
        """Getter function for magnetic moment of site"""
        return self._mag_moment

    @property
    def charge(self):
        """Getter function for charge state of site"""
        return self._charge

    @property
    def defect(self):
        """Getter function for charge state of site"""
        return self._isdefect

class Structure:
    def __init__(self, lattice, species=None,
     coords=None, sites=None, basis=None, space_group=None, defects=None, scale_matrix=None, name=None):

        """Class Structure defined by collective set of sites and lattice. Sites can either be provided with species list and coords list seperatley
         or list of sites (type Sites)
        **Args:

        lattice (ndarray): lattice describing structure in crystallagraphic system
        species (list): list of species in structure [default = None]
        coords (list): list of coordinates in structure [default = None]
        sites (list): list of sites in structure [default = None]
        basis (list): list of basis sites [default = initialized list of sites]
        space_group (int): space group number describing symmetry of structure [default = None]
        defects (list): list of defect site [default = None]
        scale_matrix (ndarray): 3x3 matrix respresenting size of supercell [default = [[1,0,0],[0,1,0],[0,0,1]]]
        name (str): name of structure [default="struct"]
        """

        self._name = name or "struct"
        assert isinstance(lattice, np.ndarray)
        assert np.shape(lattice) == (3,3)
        self._lattice = lattice
        defect_lst = []
        if sites:
            itr = 1
            for site in sites:
                assert isinstance(site, Site), "site must be of type Site!!"
                site._site_idx = itr
                itr += 1
                if site._isdefect == True:
                   defect_lst.append(site)
                   continue
                else:
                   continue
            self._sites = sites

        elif species:
            assert len(coords) == len(species), "coords and species must be the same length"
            sites = []
            for i in range(len(species)):
                coord = coords[i]
                spec = species[i]
                assert isinstance(spec, str), "all species must be type string!!"
                assert isinstance(coord, np.ndarray), "all species must be type ndarray!!"
                assert len(coord) == 3
                sites.append(Site(spec, coord))
            self._sites = sites
        else:
            print("Must provide input as either sites or both coordinates and species to generate structure!!")

        if basis:
           for site in self._basis:
               assert isinstance(site, Site), "basis must be a list of Site objects!!"
           self._basis = basis

        else:
          self._basis = self._sites

        self._defects = defects or defect_lst

        if space_group:
            assert isinstance(space_group, int) and space_group >= 1 and space_group <= 230, "space group must be an integer"
            self._space_group = space_group
        else:
            self._space_group = 0

        self._dimensions = scale_matrix or np.array([[1,0,0], [0,1,0], [0,0,1]])
        self._natoms = len(self._sites)

    #non-default constructor
    # def from_sg(self, lattice, basis, space_group):
    #     from pymatgen.core.structure import from_spacegroup
    #     ## insert propoer assertions
    #     self._lattice = lattice
    #     self._basis = basis
    #     self._space_group = space_group
    #     coords = np.array([self._basis.coord for coord in self._basis])
    #     species = [self._basis.species for species in self._basis]
    #     struct = pg.core.Structure.from_spacegroup(self._space_group, lattice, species, coords)
    #     self._dimensions  = np.array([[1,0,0], [0,1,0], [0,0,1]])

    # @classmethod
    # def make_super(self, dimensions):
    #     self._dimensions = dimensions
    #     species = [self._basis.species for species in self._basis]
    #     coords = np.array([self._basis.coord for coord in self._basis])
    #     struct = pg.core.Structure(self._lattice, species, coords)
    #     struct.make_supercell(self._dimensions)
    #     self._lattice = structure.lattice()
    #     sites = struct._sites
    #     self._sites = []
    #     for site in sites:
    #         temp_site = Site(site._species, site._coords)
    #         self._sites.append(temp_site)

    # @classmethod
    # def add_site(self, site):
    #     assert type(site) is Site, "site must be of type Site"
    #     self._sites.append(site)
