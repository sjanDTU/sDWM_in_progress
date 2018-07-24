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
    #print 'ti', ti
    xi = MannBox.U * (MannBox.T - ti)
    #print 'xi', xi
    index_xi = int(xi * MannBox.nx_r / MannBox.lx) - 1
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
#Not needed
def Interpolate_with_polargrid(MannBox):
    y = np.array(np.linspace(-MannBox.ly / 2, MannBox.ly / 2, MannBox.ny))
    z = np.array(np.linspace(-MannBox.lz / 2, MannBox.lz / 2, MannBox.nz))

    r, t = cart2pol(y, z)

    # start_time = time.time()
    if MannBox.Interp_method == 'interp2d':
        fluct_y_z = interp2d(r, t, list(MannBox.plan_of_interest), kind='linear')
    if MannBox.Interp_method == 'RBS':
        # Fastest way to interpolate
        fluct_y_z = RectBivariateSpline(r, t, MannBox.plan_of_interest,
                                        # bbox=[-MannBox.ly/2,MannBox.ly/2,-MannBox.lz/2,MannBox.lz/2],
                                        kx=5,  # number of splines: can't be superior to 5
                                        ky=5,  # same for y. Much more you have spline fast is the integration after
                                        s=1)  # positive smoothing factor
    return fluct_y_z

def Interpolate_with_polargrid_v2(MannBox):
    y = np.linspace(-MannBox.ly / 2, MannBox.ly / 2, MannBox.ny)
    z = np.linspace(-MannBox.lz / 2, MannBox.lz / 2, MannBox.nz)

    yy, zz = np.meshgrid(y,z)
    print yy
    print zz

    rr, tt = cart2pol(yy, zz)

    shape = np.shape(rr)

    r_sbs=[]
    t_sbs=[]
    val_sbs=[]
    for i in range(shape[0]):
        for k in range(shape[1]):
            r_sbs.append(rr[i,k])
            t_sbs.append(tt[i,k])
            val_sbs.append(MannBox.plan_of_interest[i,k])

    # start_time = time.time()
    if MannBox.Interp_method == 'interp2d':
        fluct_y_z = interp2d(rr, tt, MannBox.plan_of_interest, kind='cubic')
    if MannBox.Interp_method == 'SBS':
        # Fastest way to interpolate
        fluct_y_z = SmoothBivariateSpline(r_sbs, t_sbs, val_sbs,
                                        # bbox=[-MannBox.ly/2,MannBox.ly/2,-MannBox.lz/2,MannBox.lz/2],
                                        kx=3,  # number of splines: can't be superior to 5
                                        ky=3,  # same for y. Much more you have spline fast is the integration after
                                        s=10)  # positive smoothing factor
        plt.figure()
        plt.subplot(121)
        plt.pcolor(yy,zz,MannBox.plan_of_interest)
        plt.subplot(122)
        plt.pcolor(yy,zz,fluct_y_z(rr,tt))
        plt.show()

    return fluct_y_z

def Trapz_for_Integrate_for_regulargrid(Interpo_Integrate):
    r = Interpo_Integrate.r
    t = Interpo_Integrate.t
    value = Interpo_Integrate.ValMesh

    Final_Integral_Value = np.trapz([np.trapz(value[i, :], r) for i in range(Interpo_Integrate.nt)], t)
    return Final_Integral_Value
def Trapz_for_Integrate_for_regulargrid(Interpo_Integrate):
    r = Interpo_Integrate.r
    t = Interpo_Integrate.t
    value = Interpo_Integrate.ValMesh

    Final_Integral_Value = np.trapz([np.trapz(value[i, :], r) for i in range(Interpo_Integrate.nt)], t)
    return Final_Integral_Value

def PlanValues_multiplied_by_radius(MannBox):
    y_vec = np.linspace(-MannBox.ly / 2, MannBox.ly / 2, MannBox.ny)
    z_vec = np.linspace(-MannBox.lz / 2, MannBox.lz / 2, MannBox.nz)

    yy, zz = np.meshgrid(y_vec, z_vec)

    #rr = cart2pol(yy, zz)[0]
    rr = np.sqrt(yy**2+zz**2)

    #MannBox.radius_plan_of_interest = np.zeros(np.shape(MannBox.plan_of_interest))

    MannBox.radius_plan_of_interest = rr * MannBox.plan_of_interest
    return MannBox
def BoxValues_multiplied_by_radius(MannBox):
    y_vec = np.linspace(-MannBox.ly / 2, MannBox.ly / 2, MannBox.ny)
    z_vec = np.linspace(-MannBox.lz / 2, MannBox.lz / 2, MannBox.nz)

    yy, zz = np.meshgrid(y_vec, z_vec)

    # rr = cart2pol(yy, zz)[0]
    rr = np.sqrt(yy ** 2 + zz ** 2)

    for i in range(MannBox.nx_r):
        MannBox.v_TurbBox[:, :, i] = rr * MannBox.v_TurbBox[:, :, i]
        MannBox.w_TurbBox[:, :, i] = rr * MannBox.w_TurbBox[:, :, i]
    return MannBox

def Integrate(MannBox, Interpo_Integrate):
    Rw = MannBox.WakeRadius_for_this_Plan

    nr = Interpo_Integrate.nr
    nt = Interpo_Integrate.nt
    r = np.linspace(0., 2 * Rw, nr)  # centre (yb,zb) and with a radius equal to Dw. (2*Rw)
    t = np.linspace(-m.pi, m.pi, nt)

    PolMesh = np.meshgrid(r, t)  # rr, tt = meshgrid(r, t)
    CartMesh = pol2cart_for_mesh(PolMesh)  # yy, zz from rr, tt

    ###### more time expansive part of the function
    if Interpo_Integrate.Interp_method == 'RBS':
        ValMesh = PolMesh[0] * Interpo_Integrate.f_cart(CartMesh[0], CartMesh[1])

        # Interpo_Integrate.ValMesh = Interpo_Integrate.f_cart(Interpo_Integrate.CartMesh[0], Interpo_Integrate.CartMesh[1])

    if Interpo_Integrate.Interp_method == 'interp2d':
        Interpo_Integrate.ValMesh = np.zeros(np.shape(Interpo_Integrate.CartMesh[0]))  # matrix of r*f_cart(y,z)
        # """
        for i in range(np.shape(Interpo_Integrate.ValMesh)[0]):
            for k in range(np.shape(Interpo_Integrate.ValMesh)[1]):
                y = Interpo_Integrate.CartMesh[0][i, k]
                z = Interpo_Integrate.CartMesh[1][i, k]

                # Calculate r*f(x,z) to be ready to be integrated
                Interpo_Integrate.ValMesh[i, k] = np.sqrt((y) ** 2 + (z) ** 2) * Interpo_Integrate.f_cart(y, z)
                # """

    I = np.zeros(nt)
    for i in range(nt):
        I[i] = np.trapz(ValMesh[i, :], r)
    return np.trapz(I, t)



########################################################################################################################

def Filter_Turb_Box(Rw, MannBox):
    """
    Function to get the two component v_c and w_c characteristic of the meandering phenomena.
    The filter depends on the wake radius, so we must know the wake radius before filtering

    (Larsen, G. C. 2008 - Wake meandering - a pragmatic approach)

    For a question of memory cost (to store the Turb_box), this code deal with one by one component.
    So you have to run it twice to get the two components (vc, wc)
    :param Turb_box:
    :return:
    """
    MannBox.plan_of_interest = get_plane_of_interest(MannBox, ti=0.)

    f = Interpolate_plan(MannBox)

    """ Initialize the algo """

    # Wake center location component related to the turb_Box concerned (ti=0, x=0)
    # (y_g for v_c and so for v_fluct and z_g for w_c and so for w_fluct)
    ti = 0.  # (m/s)
    yi = 0.;
    zi = 0.  # (m)
    #Af_index = get_Af_for_filtering(Rw, ti, yi, zi, MannBox)

    xi = MannBox.U * (MannBox.T - ti)
    index_xi = int(xi * MannBox.nx / MannBox.lx) - 1

    Plan_of_interest = MannBox.TurbBox[:, :, index_xi]
    MannBox.TurbBox = MannBox.TurbBox[:, :, :-1]  # delete one plan in the turbbox to
    # reduce memory cost along the algorithm

    #DATA = extract_data_for_integ_2d(Plan_of_interest, Af_index)

    # Determine v_c, w_c characteristic velocities from Turb


    # Define the averaged area:
    # The selected averaging area, Af, is most logically selected as a circle concentric
    # with the advected wake deficit (i.e. with centre (yb,zb) and with a radius equal to Dw.
    Af = m.pi * Rw ** 2.

    # Loop for double Integration
    Af = m.pi * Rw ** 2.

    return

########################################################################################################################

#print dblquad(lambda t, r: r,0,1, lambda t: 0, lambda t: 2*m.pi)

#func = lambda r,t : r
#options={'limit':20,'epsrel':10**-6, 'epsabs':10**-6}
#integral=nquad(func,[[0,1],[0,2*m.pi]],opts=[options,options])
#print integral