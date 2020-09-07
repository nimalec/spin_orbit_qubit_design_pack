import numpy as np
from src.task import *

class BWOFe_MAEworkflow(MagenticAnisotropySphereFlow):

    def __init__(self, workdir, npoints, nbands, struct_path, ref_orient=[0,0,1], cl_path=None, nelect=None, name="BWOFe_mae_sphere"):
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

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path, name ="BWOFe_mae_sphere")



class PTOFe_MAEworkflow(MagenticAnisotropySphereFlow):

    def __init__(self, workdir, npoints, nbands, struct_path, ref_orient=[0,0,1], cl_path=None, nelect=None, name="PTOFe_mae_sphere"):
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
<<<<<<< HEAD:example_setup_scripts/example_complexoxide.py
            potcar_path_ = "../pseudos/PTO_Fe_POTCAR"
=======
            potcar_path_ = "/global/scratch/nleclerc/spin_orbit_qubit_design_pack/pseudos/PTOFe_POTCAR" 
>>>>>>> 336b2815fa1b34a17b5a04ad8a8fcea422e6efe9:src/example_complexoxide.py
            kgrid_ = [2,2,2]
            nodes_ = 6
            ppn_ = 24
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0, 0, 0, 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 750
            magmom_ = [0, 0, 0, 6]

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path, name ="PTOFe_mae_sphere", ismear=0)
