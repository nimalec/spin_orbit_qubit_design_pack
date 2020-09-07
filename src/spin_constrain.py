import numpy as np
from task import *

class BWOFeSpinConstr(SpinMultipletSeries):

    def __init__(self, workdir, npoints, nupdown_rng,struct_path, ref_orient=[0,0,1],nbands=None, nelect=None, name="BWOFe_spin"):
            """
            Computes MCAE for Fe doped Bi2WO6 materials.

            **Args:

            workdir (str):
            npoints (int):
            pseudo_dir (str):
            struct_path (str):
            cl_path (str):

            """
            self._workdir = workdir
            self._npoints = npoints
            potcar_path_ = "../pseudos/BWO_Fe_POTCAR"
            kgrid_ = [2,2,2]
            nodes_ = 6
            ppn_ = 24
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0, 0, 0, 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 800
            magmom_ = [0, 0, 0, 6]


            SpinMultipletSeries.__init__(self, workdir=workdir, npoints=npoints, nupdown_rng=nupdown_rng ,kgrid=kgrid_, nodes=nodes_, ppn=ppn_, ref_orient=ref_orient, ldaul=ldaul_, magmom=magmom_, Uparam=Uparam_, Jparam=Jparam_, encut=encut_, potcar_path=potcar_path_, struct_path=struct_path, ismear=-5, sigma=0.2, nelect=nelect, nbands=nbands)
