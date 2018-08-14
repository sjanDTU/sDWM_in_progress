"""

"""
import numpy as np
from scipy.interpolate import interp2d, RectBivariateSpline
from scipy.integrate import trapz
import math as m
import matplotlib.pyplot as plt
from Plot_Function import *
from cTurb import interpo_integrate
import time

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y
def pol2cart_for_mesh(mesh):
    x = mesh[0] * np.cos(mesh[1])
    y = mesh[0] * np.sin(mesh[1])
    return x, y
def cart2pol_for_mesh(mesh):
    rho = np.sqrt(mesh[0]**2 + mesh[1]**2)
    phi = np.arctan2(mesh[1], mesh[0])
    return rho, phi

def get_plan_of_interest(MannBox,ti,component_char):
    # print 'ti', ti
    # print MannBox.U
    # print MannBox.T
    # print MannBox.lx
    xi = MannBox.U * (MannBox.T - ti)
    #print 'xi', xi
    if MannBox.nx_r!=None:
        index_xi = int(xi * MannBox.nx_r / MannBox.lx) - 1
    else:
        index_xi = int(xi * MannBox.nx / MannBox.lx) - 1
    #print 'index_xi: ', index_xi
    if component_char == 'vfluct':
        return MannBox.v_TurbBox[:, :, index_xi]

    if component_char == 'wfluct':
        return MannBox.w_TurbBox[:, :, index_xi]
        #MannBox.w_TurbBox = MannBox.w_TurbBox[:, :, :-1]


def Interpolate_plan(MannBox, Interpo_Integrate, Meand_Mann):

    yc, zc = Meand_Mann.wake_center_location

    # Define a referential Wakecentered
    y = np.linspace(-MannBox.ly/2-yc, MannBox.ly/2-yc, MannBox.ny)
    z = np.linspace(-MannBox.lz/2-zc, MannBox.lz/2-zc, MannBox.nz)


    #start_time = time.time()
    if Interpo_Integrate.Interp_method == 'interp2d':
        Interpo_Integrate.f_cart = interp2d(y, z, MannBox.plan_of_interest, kind = Interpo_Integrate.Interp2dKind, copy = True)
    if Interpo_Integrate.Interp_method == 'RBS':
        # Fastest way to interpolate
        Interpo_Integrate.f_cart = RectBivariateSpline(x=y, y=z, z=np.transpose(MannBox.plan_of_interest),  # np.transpose for last version of scipy
                                    kx=1,   # number of splines: can't be superior to 5
                                    ky=1,
                                    s=0)    # positive smoothing factor
    return Interpo_Integrate


def polar_mesh_value(MannBox, Interpo_Integrate):

    Rw = MannBox.WakeRadius_for_this_Plan
    #Rw = MannBox.ly/4 #To plot a check with a larger aera

    nr = Interpo_Integrate.nr
    nt = Interpo_Integrate.nt
    Interpo_Integrate.r = np.linspace(0., 2*Rw, nr)   # centre (yb,zb) and with a radius equal to Dw. (2*Rw)
    Interpo_Integrate.t = np.linspace(-m.pi, m.pi, nt)

    Interpo_Integrate.PolMesh = np.meshgrid(Interpo_Integrate.r, Interpo_Integrate.t)                              # rr, tt = meshgrid(r, t)
    Interpo_Integrate.CartMesh = pol2cart_for_mesh(Interpo_Integrate.PolMesh)  #yy, zz from rr, tt

    ###### more time expansive part of the function
    if Interpo_Integrate.Interp_method == 'RBS':

        # The following line need a change in the fit2pack package: L774 in def __call__ put ---> grid=False
        Interpo_Integrate.ValMesh = Interpo_Integrate.PolMesh[0] * Interpo_Integrate.f_cart(Interpo_Integrate.CartMesh[0], Interpo_Integrate.CartMesh[1])

    elif Interpo_Integrate.Interp_method == 'interp2d':
        Interpo_Integrate.ValMesh = np.zeros(np.shape(Interpo_Integrate.CartMesh[0]))  # matrix of r*f_cart(y,z)
        #"""
        for i in range(np.shape(Interpo_Integrate.ValMesh)[0]):
            for k in range(np.shape(Interpo_Integrate.ValMesh)[1]):
                y = Interpo_Integrate.CartMesh[0][i, k]
                z = Interpo_Integrate.CartMesh[1][i, k]

                # Calculate r*f(x,z) to be ready to be integrated
                Interpo_Integrate.ValMesh[i, k] = np.sqrt((y)**2+(z)**2)*Interpo_Integrate.f_cart(y, z)
    else:
        print 'Interp Method is not recognized'

    if Interpo_Integrate.Plot:
        plot_interpolation(MannBox, Interpo_Integrate)

    return Interpo_Integrate


def Trapz_for_Integrate_general_grid(Interpo_Integrate):
    r = Interpo_Integrate.r
    t = Interpo_Integrate.t
    value = Interpo_Integrate.ValMesh

    Final_Integral_Value = trapz([trapz(value[i, :], r) for i in range(Interpo_Integrate.nt)], t)
    #Final_Integral_Value = np.trapz([trapz(value[i, :], r) for i in range(Interpo_Integrate.nt)], t)

    return Final_Integral_Value




########################################################################################################################
#### if we need a complex interp
# Work in progress if needed!
# but we can just do a simple np.interp
def interpo_vc_wc_time(MannBox,Meand_Mann, Interpo_Integrate):
    print 'shape init vc wc: ', np.shape(Meand_Mann.init_vc_wc)
    vc = Meand_Mann.init_vc_wc[:,0]
    wc = Meand_Mann.init_vc_wc[:,1]

    T_disc = np.array(MannBox.ti)
    print 'shape t_disc: ', np.shape(T_disc)

    Interpo_Integrate.f_vc_t = np.interp()
    return





########################################################################################################################


########################################################################################################################



########################################################################################################################

#print dblquad(lambda t, r: r,0,1, lambda t: 0, lambda t: 2*m.pi)

#func = lambda r,t : r
#options={'limit':20,'epsrel':10**-6, 'epsabs':10**-6}
#integral=nquad(func,[[0,1],[0,2*m.pi]],opts=[options,options])
#print integral