import numpy as np
from task import *

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

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path, ismear='0', name ="BWOFe_mae_sphere", nelect=nelect, sigma=0.05)

class BWOFe_W_MAEworkflow(MagenticAnisotropySphereFlow):

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
            sigma_ = 0.05    
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0, 0, 0, 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 800
            magmom_ = [0, 0, 0, 6]

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path, ismear=str(0), name ="BWOFe_mae_sphere", nelect=nelect,  sigma=sigma_)  


class BWOCo_MAEworkflow(MagenticAnisotropySphereFlow):

    def __init__(self, workdir, npoints, nbands, struct_path, ref_orient=[0,0,1], cl_path=None, nelect=None, name="BWOCo_mae_sphere"):
            """
            Computes MCAE for Mn doped Bi2WO6 materials.

            **Args:

            workdir (str):
            npoints (int):
            pseudo_dir (str):
            struct_path (str):
            cl_path (str):

            """
            self._workdir = workdir
            self._npoints = npoints
            potcar_path_ = "../pseudos/BWO_Co_POTCAR"
            kgrid_ = [2,2,2]
            nodes_ = 6
            ppn_ = 24
            sigma_ = 0.05    
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0, 0, 0, 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 800
            magmom_ = [0, 0, 0, 6]

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path, ismear=str(0), name ="BWOCo_mae_sphere", nelect=nelect,  sigma=sigma_)  


class BWOMn_MAEworkflow(MagenticAnisotropySphereFlow):

    def __init__(self, workdir, npoints, nbands, struct_path, ref_orient=[0,0,1], cl_path=None, nelect=None, name="BWOMn_mae_sphere"):
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
            potcar_path_ = "../pseudos/BWO_Mn_POTCAR"
            kgrid_ = [2,2,2]
            nodes_ = 6
            ppn_ = 24
            sigma_ = 0.05    
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0, 0, 0, 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 800
            magmom_ = [0, 0, 0, 6]

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path, ismear=str(0), name ="BWOMn_mae_sphere", nelect=nelect,  sigma=sigma_)  

class PTOFe_MAEworkflow(MagenticAnisotropySphereFlow):

    def __init__(self, workdir, npoints, nbands, struct_path, ref_orient=[0,0,1], ang_rng=[[0,0.5*np.pi],[0, np.pi]], cl_path=None, nelect=None, name="PTOFe_mae_sphere"):
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
            potcar_path_ = "../pseudos/PTOFe_POTCAR"
            kgrid_ = [2,2,2]
            nodes_ = 6
            ppn_ = 24
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0, 0, 0, 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 750
            magmom_ = [0, 0, 0, 6]

            MagenticAnisotropySphereFlow.__init__(self, workdir, npoints, kgrid_, nbands, nodes_, ppn_, ref_orient, ldaul_, magmom_, Uparam_, Jparam_, encut_, potcar_path_, struct_path,angle_rng=ang_rng,ismear='0', name ="BWOFe_mae_sphere", nelect=nelect, sigma=0.01)
