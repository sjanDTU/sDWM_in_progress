"""
date: 24/07/18
author: Augustin Grosdidier
Purpose: Treatment of a MannBox for wake added turbulence process
"""

import numpy as np
import math
import matplotlib.pyplot as plt

from Polar_Interpolate_Data_Plan import cart2pol_for_mesh
from ReadMannTurbulence import *

from scipy.interpolate import griddata, RectBivariateSpline



def pre_init_turb_WaT(filename):
    """
    Purpose:
    Extract Ti_u for wake added Components
    :return:
    """

    # ------------------ # GET TURBULENT BOX # -------------------------- #
    MannBox = ReadMannInput(filename)

    MannBox.U_ref = get_averaged_U(filename)  # Specified U_mean for this Mannbox

    MannBox.u_TurbBox = get_turb_component_from_MannBox(filename, 'ufluct', False, MannBox, video=False)/MannBox.U_ref
    MannBox.TI_u = np.sqrt(np.mean((MannBox.u_TurbBox) ** 2))
    print 'TI from MannBox WaT: ', MannBox.TI_u

    # Size the MannBox according to MFOR / FFoR / mixl domain settings

    return MannBox

def init_turb_WaT(MannBox, meta):
    """

    run at each loop iteration, scale mannbox onto mean turbulence intensity

    :param MannBox:
    :return:
    """

    # is it possible to not create a second mannbox to save memory cost
    # ---------------------------- # Scale Turbulence Intensity # -------------------- #
    # We take of the last scaling due to the last iteration
    MannBox.u_TurbBox = MannBox.u_TurbBox / MannBox.k_scale

    # we update/calculate the new scaling for the next iteration
    MannBox.k_scale = meta.mean_TI_DWM / MannBox.TI_u

    # We apply the scaling
    MannBox.u_TurbBox = MannBox.k_scale * MannBox.u_TurbBox

    return

def obtain_wake_added_turbulence(MannBox, i_t, meta):
    """
    Get the plan of wake added Turbulence
    Prepare it to be expressed in the meandering Frame
    Apply km_t(r)
    return the wake added turbulence plan ready to be added in the meandering frame
    Notice:
    it will be expressed in FFoR as the deficit coming from Ainslie.

    :param MannBox:
    :param meta:
    :return:
    """

    plan_of_interest = MannBox.u_TurbBox[:, :, i_t]
    Ly = np.linspace(meta.x_vec[0], meta.x_vec[-1], MannBox.ny)
    Lz = np.linspace(meta.y_vec[0], meta.y_vec[-1], MannBox.nz)
    # we consider this plan centered on the wake
    #new_ly = np.linspace(meta.x_vec[0], meta.x_vec[-1], meta.nx)
    #new_lz = np.linspace(meta.y_vec[0], meta.y_vec[-1], meta.ny)
    new_ly = meta.x_vec
    new_lz = meta.y_vec

    #WaT = griddata((Ly,Lz),plan_of_interest,(new_ly, new_lz))
    #WaT = RectBivariateSpline()

    f_cart = RectBivariateSpline(x=Ly, y=Lz, z=np.transpose(plan_of_interest),
                                                   # np.transpose for last version of scipy
                                                   kx=1,  # number of splines: can't be superior to 5
                                                   ky=1,
                                                   s=0)  # positive smoothing factor

    yy, zz = np.meshgrid(new_ly, new_lz)

    WaT = f_cart(yy,zz)

    # reshape WaT to be compatible with FFoR
    """
    plt.figure()
    plt.contourf(Distrib*plan_of_interest)
    plt.colorbar()
    plt.show()
    #"""
    return WaT





"""
# Test
#result = pre_init_turb_WaT('1028')
#print result.TI_u

kmt_r = np.load('kmt_r_for_iz0.NPY')
vr_mixl =  np.load('vr_mixl.NPY')
plt.figure()
plt.plot(vr_mixl,kmt_r)
plt.show()
print kmt_r
print vr_mixl

MannBox = pre_init_turb_WaT('1101')
obtain_wake_added_turbulence(MannBox, 50, vr_mixl, kmt_r)
#"""