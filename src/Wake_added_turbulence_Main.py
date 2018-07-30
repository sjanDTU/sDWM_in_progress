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
    MannBox.U_ref = get_averaged_U(filename)# Specified U_mean for this Mannbox
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

def obtain_wake_added_turbulence(MannBox, i_t, vr_mixl, kmt_r):
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
    # we consider this plan centered on the wake

    """
    plt.figure()
    plt.contourf(plan_of_interest)
    plt.colorbar()
    plt.show()
    #"""

    #yy, zz = np.meshgrid(np.linspace(-meta.vr_mixl, meta.vr_mixl, MannBox.ny),
     #                    np.linspace(-meta.vr_mixl, meta.vr_mixl, MannBox.nz))
    yy, zz = np.meshgrid(np.linspace(-vr_mixl[-1], vr_mixl[-1], MannBox.ny),
                        np.linspace(-vr_mixl[-1], vr_mixl[-1], MannBox.nz))
    rr = cart2pol_for_mesh((yy, zz))[0]

    Distrib = np.interp(rr, vr_mixl, kmt_r)
    WaT = Distrib * plan_of_interest

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