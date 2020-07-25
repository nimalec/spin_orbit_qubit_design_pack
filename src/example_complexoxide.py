import numpy as np
from task import *

class BWOFe_MAEworkflow(MagenticAnisotropySphereFlow):

    def __init__(self, workdir, npoints, nbands, potcar_path, struct_path, ref_orient, cl_path=None, name="BWO_calc"):
            """
            Computes MCAE for Fe doped Bi2WO6 materials.

            **Args:

            workdir (str):
            npoints (int):
            pseudo_dir (str):
            struct_path (str):
            cl_path (str):

            """
            self._workdir = work_dir
            self._npoints = npoints
            kgrid_ = [2,2,2]
            nodes_ = 6
            ppn_ = 24
            ldaul_ =  [-1, -1, -1, 2]
            Uparam_ = [0 0 0 4]
            Jparam_ = [0, 0, 0, 0]
            encut_ = 800
            magmom_ = [0, 0, 0, 6]

            MagenticAnisotropySphereFlow.__init__(self, workdir=self._workdir, npoints=self._npoints, kgrid=kgrid_, nbands=nbands, nodes=nodes_, ppn=ppn_, ref_orient=ref_orient, ldaul=ldaul_, name=name, magmom=magmom_,
            Uparam=Uparam_, Jparam=Jparam_, encut=encut_, potcar_path=potcar_path, struct_path=struct_path)

if __name__ = "main":
    BWO_Fe_Bi_poscar_path_ = "/global/scratch/nleclerc/soc_mae_predictions/bwo_fe_maesphere_v0/Fe_Bi"
    #BWO_Fe_W_workdir_ =
    npoints_ = 200
    nbands_ = 1494
    potcar_path_ = "/global/scratch/nleclerc/soc_mae_predictions/pseudos/BiWOFe/BiWFeO_POTCAR"
    BWO_Fe_Bi_poscar_path_ = "/global/scratch/nleclerc/soc_mae_predictions/test_mae_v3/test_cl/POSCAR"
    #BWO_Fe_W_poscar_path_ =
    ref_orient_ = [0,0,1]
    mae_calc_BWO_Fe_Bi_ = BWOFe_MAEworkflow(BWO_Fe_W_workdir_, npoints_, nbands, potcar_path, BWO_Fe_Bi_poscar_path_, ref_orient_)
    #mae_calc_BWO_Fe_W_  = BWOFe_MAEworkflow(BWO_Fe_W_workdir_, npoints_, nbands, potcar_path, BWO_Fe_Bi_poscar_path_, ref_orient)
    mae_calc_BWO_Fe_Bi_.make_calculations()
    #mae_calc_BWO_Fe_W_.make_calculations()
