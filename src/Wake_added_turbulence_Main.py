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
    print 'MannBoxsize: ('+str(MannBox.ly)+', '+str(MannBox.lz)+')'
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
    debug = True
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
    MannBox.k_scale = meta.mean_TI_DWM / (MannBox.TI_u * MannBox.U_ref / meta.WS)

    # We apply the scaling
    MannBox.u_TurbBox = MannBox.k_scale * MannBox.u_TurbBox
    if debug:
        print 'kscale: ', MannBox.k_scale
        print 'Mannbox Ti', np.sqrt(np.mean((MannBox.u_TurbBox) ** 2))/ meta.WS
        #raw_input('...')

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
    if meta.Meandering:
        t = meta.time[i_t] # needed to synchronize in time Meandering and WaT
        xi =  meta.WS / MannBox.U  * t  # get the plan of interest , scale to the correct ambien WS
        index_t = int(round(xi))

    plan_of_interest = MannBox.u_TurbBox[:, :, index_t]

    # Box containing cross section covering one rotor diameter (according to Madsen, 2010, calibration ...)
    if MannBox.One_rotordiameter_size:
        Ly = np.linspace(-MannBox.R_MB, MannBox.R_MB, MannBox.ny)
        Lz = np.linspace(-MannBox.R_MB, MannBox.R_MB, MannBox.nz)
    elif MannBox.based_on_MannBoxsize:
        Ly = np.linspace(-MannBox.ly/2, MannBox.ly/2, MannBox.ny)
        Lz = np.linspace(-MannBox.lz / 2, MannBox.lz / 2, MannBox.ny)
        Ly = Ly /MannBox.R_ref
        Lz = Lz/MannBox.R_ref

    f_cart = RectBivariateSpline(x=Ly, y=Lz, z=np.transpose(plan_of_interest/meta.WS),
                                                   # np.transpose for last version of scipy
                                                   kx=1,  # number of splines: can't be superior to 5
                                                   ky=1,
                                                   s=0)  # positive smoothing factor

    return f_cart






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