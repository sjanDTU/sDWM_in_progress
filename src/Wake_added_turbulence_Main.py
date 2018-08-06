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

    MannBox.u_TurbBox = get_turb_component_from_MannBox(filename, 'ufluct', False, MannBox, video=False)#/MannBox.U_ref
    MannBox.TI_u = np.sqrt(np.mean((MannBox.u_TurbBox) ** 2)) /MannBox.U_ref #RMS
    print 'TI from MannBox WaT: ', MannBox.TI_u
    print 'RMS TI, at a location: ', np.mean(np.sqrt(np.mean((MannBox.u_TurbBox) ** 2, axis = 2)) )/MannBox.U_ref #Have to be similar to MannBOx Ti_u to be Homogeneous
    """
    #Variance
    Var = np.mean((MannBox.u_TurbBox-np.mean(MannBox.u_TurbBox))**2)
    print Var
    raw_input('...')
    #"""

    # Size the MannBox according to MFOR / FFoR / mixl domain settings

    return MannBox

def init_turb_WaT(MannBox, meta):
    """

    run at each loop iteration, scale mannbox onto mean turbulence intensity

    :param MannBox:
    :return:
    """
    debug = False
    # is it possible to not create a second mannbox to save memory cost
    # ---------------------------- # Scale Turbulence Intensity # -------------------- #
    # We take of the last scaling due to the last iteration
    MannBox.u_TurbBox = MannBox.u_TurbBox / MannBox.k_scale

    # we update/calculate the new scaling for the next iteration
    if debug:
        print 'First scaling: '

        print 'k_scale', MannBox.k_scale
        print 'MannBox_Tiu', MannBox.TI_u
        print 'mean DWM Ti', meta.mean_TI_DWM
    MannBox.k_scale = meta.mean_TI_DWM / MannBox.TI_u

    # We apply the scaling
    MannBox.u_TurbBox = MannBox.k_scale * MannBox.u_TurbBox
    if debug:
        print 'kscale: ', MannBox.k_scale
        print 'Mannbox Ti', np.sqrt(np.mean((MannBox.u_TurbBox) ** 2))
        raw_input('...')

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

    # Box containing plan of (2R, 2R)
    Ly = np.linspace(0., 4., MannBox.ny)
    Lz = np.linspace(0, 4., MannBox.nz)


    new_ly = meta.x_vec     #more than 2
    new_lz = meta.y_vec


    f_cart = RectBivariateSpline(x=Ly, y=Lz, z=np.transpose(plan_of_interest),
                                                   # np.transpose for last version of scipy
                                                   kx=1,  # number of splines: can't be superior to 5
                                                   ky=1,
                                                   s=0)  # positive smoothing factor

    yy, zz = np.meshgrid(new_ly, new_lz)

    WaT = f_cart(yy,zz)/MannBox.U_ref

    #plt.figure()
    #plt.contourf(WaT)
    #plt.colorbar()
    #plt.show()

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
result = pre_init_turb_WaT('1101')
print result.TI_u
raw_input('...')


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