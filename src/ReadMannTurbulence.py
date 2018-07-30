"""
LOAD MANN BOXES
author: Augustin Grosdidier
date: 31/05/2018 at 15:01
"""
import numpy as np
import matplotlib.pyplot as plt
import math as m
from cMann import MannBox

from matplotlib import animation



def ReadMannInput(filename):
    # Mann class object Init
    Mannbox = MannBox()

    # Open datafile
    dir = 'C:/Users/augus/Documents/Stage/ClusterDATA/Mann/'
    src = dir + filename + '/mann.inp'
    f = open(src)
    Input = f.read()
    new_Input=[]
    loop_input=''
    for char in Input:
        if char!='\n':
            loop_input = loop_input + char

        else:
            new_Input = new_Input + [loop_input]
            loop_input = ''
    Input = new_Input

    print 'MannBox read for Wake added Turubulence'
    print Input

    Mannbox.fieldDim = int(Input[0]); Mannbox.N_Comp = int(Input[1])

    # Cas 'Basic'


    Mannbox.nx = int(Input[2]); Mannbox.ny = int(Input[3]); Mannbox.nz = int(Input[4])
    Mannbox.lx = float(Input[5]); Mannbox.ly = float(Input[6]); Mannbox.lz = float(Input[7])


    Mannbox.comSpec = Input[8]
    Mannbox.alpha_epsilon_2_3_spec = float(Input[9])
    Mannbox.L_spec = float(Input[10])
    Mannbox.Gamma_spec = float(Input[11])

    Mannbox.U = get_averaged_U(filename)
    Mannbox.L = Mannbox.lx
    Mannbox.T = Mannbox.L / Mannbox.U

    Mannbox.dx = Mannbox.lx / Mannbox.nx
    Mannbox.dt = Mannbox.dx / Mannbox.U

    Mannbox.ti = np.linspace(0., Mannbox.T, Mannbox.nx)

    return Mannbox

def get_averaged_U(filename):
    """
    Average Fluct NOP

    Autre Piste a partir de Mann.inp,
    {\alpha*\epsilon^{2/3}, L, \Gamma} = {0.0605, 26.1, 2.99} , <U>=11.7 m/s
    :param filename:
    :param Input:
    :return:
    """
    if filename=='1028':
        # {\alpha*\epsilon^{2/3}, L, \Gamma} = {0.0605, 26.1, 2.99} , <U>=11.7 m/s , <dir>=9.8 degree ;
        return 11.7
    if filename=='1101':
        # {\alpha*\epsilon^{2/3}, L, \Gamma} = {0.0553, 55.7, 3.81}, <U>=10.7 m/s , <dir>=13.6 degree ;
        return 10.7
    return None

############################ ADAPT IT for Wake added Turbulence
def sizing_MannBox(MannBox , WindFarm):
    """
    Note: this sizing imply a choice of the reference wich is obvious if all turbines are same => there is ONE Radius
    :param MannBox:
    :param R_WindTurbine:
    :param U_mean_WindFarm:
    :param windFarm_lenght:
    :return:
    """

    MannBox.U = MannBox.U / WindFarm.U_mean
    print 'U (no dimension): ', MannBox.U
    MannBox.L = MannBox.L / WindFarm.WT_R * k
    print 'L (no dimension): ', MannBox.L

    MannBox.T = MannBox.L / MannBox.U; print 'Mannbox T_total (no dimension): ', MannBox.T
    MannBox.dt = MannBox.T / MannBox.nx; print 'dt (no dimension): ', MannBox.dt
    MannBox.ti = [i*MannBox.dt for i in range(MannBox.nx)]
    MannBox.ti = [t for t in MannBox.ti
                  if t < MannBox.SimulationTime + MannBox.dt]  # we catch just one point after the simulation time


    MannBox.lx = MannBox.L
    MannBox.ly = MannBox.ly / WindFarm.WT_R * k
    MannBox.lz = MannBox.lz / WindFarm.WT_R * k

    MannBox.dx = MannBox.lx / MannBox.nx
    MannBox.dt = MannBox.dx / MannBox.U

    MannBox.U_ref = WindFarm.U_mean
    MannBox.R_ref = WindFarm.WT_R

    MannBox.u_TurbBox = MannBox.u_TurbBox / MannBox.U_ref
    return MannBox

def get_turb_component_from_MannBox(filename,kind_of_fluct,plot_bool,MannBox, video):
    """
    It is originally a Matlab code to read Mann Turbulence.
    ReadMannTurbulence.m created by Soren Andersen (May 14th 2012)
    Purpose:
    Extract the fluctuating speed (u', v', w') in the plane (x-U*t, y, z, 0) meaning along x-axis
    from the total wind speed field (u', v', w') in (x, y, z, t) referential.
    We just extract one of the component to avoid to crowd computer ressources.
    Therefore to extract all the component, you need to use this function then a filter on the component extracted,
    then repeat it for the two other components

    :param filename:
    :param kind_of_fluct: ufluct (u-component), vfluct (v-component), wfluct (w-component)
    notice: u-component is not used for DWM.
    :param plot_bool:
    :return:
    """
    plt.close('all')

    dir = 'C:/Users/augus/Documents/Stage/ClusterDATA/Mann/'
    src = dir + filename+'/'

    """ Find these values in the mann.inp """
    nx = MannBox.nx; ny = MannBox.ny; nz = MannBox.nz         # Number of points in each direction
    lx = MannBox.lx; ly = MannBox.ly; lz = MannBox.lz         # Length in each direction[m]

    # For non dimension-result
    #U0 = 10.                                 # Freestream velocity
    #R = 40.                                  # Turbine radius
    """ Loading Mann boxes """
    count = nx*ny*nz

    fid = open(src+kind_of_fluct)
    Total_Vel = np.fromfile(fid, dtype=np.float32, count=-1)


    """ Reshaping to 3D box"""
    Total_Vel = np.reshape(Total_Vel,(ny,nz,nx), order='F')
    Total_Vel = Total_Vel[:, :, ::MannBox.discr_reduc_factor]
    # /!\ with dimension /!\
    # (it will turn to no dimension figures with sizing function wich depends of TI)
    # that's why we don't turn to no dimension right now

    # MannBox.discr_reduc_factor: discretization reducing, if you want to just compute for each 2,3,...,k plans
    # Otherwise it's 1

    """ Y and Z vectors """
    Y = np.linspace(-ly/2., ly/2., ny)
    Z = np.linspace(-lz/2., lz/2., nz)

    """ Plot Part """
    #"""
    if plot_bool:
        #fluc = np.empty((ny, nz, nx))
        #fluc[:] = np.nan

        plt.close('all')
        plt.ion()
        plt.figure()

        for i in range(nx):

            Max_fluc = np.max(Total_Vel[:, :, i])
            Min_fluc = np.min(Total_Vel[:, :, i])
            levels=np.linspace(Min_fluc, Max_fluc, 6)
            plt.clf()
            plt.cla()
            plt.contourf(Y, Z, Total_Vel[:, :, i], levels=levels, cmap = plt.cm.jet)
            plt.colorbar(), plt.xlabel('y [R]'), plt.ylabel('z [R]'),plt.title(kind_of_fluct)
            plt.draw()
            plt.pause(0.1)
    #plt.show()
    #plt.ioff()
    #"""
    if video:
        T_video = 10
        number_of_frames = int(T_video/MannBox.dt)  # for 10s: 200 frames
        print 'number_of frame'
        # Figure Definition
        fig = plt.figure()
        fig.suptitle('Mannbox for '+kind_of_fluct)
        plt.xlabel('y'), plt.ylabel('z')

        # Image Definition

        ims = []
        for i in range(number_of_frames):
            im = plt.imshow(Total_Vel[:, :, i],
                            #cmap=plt.cm.jet,
                            extent=(-MannBox.ly/2, MannBox.ly/2, -MannBox.lz/2, MannBox.lz/2),
                            animated = True)  # voir extent pour axe
            ims.append([im])
        ani = animation.ArtistAnimation(fig, ims, interval=MannBox.dt*1000)

        ani.save('C:/Users/augus/Documents/Stage/Presentation/Video/Mannbox.html')
        plt.show()
        print 'Video process ended'

    print 'shape(fluc_component_field): ', np.shape(Total_Vel)
    return Total_Vel

def sizing_MannBox_MFOR(MannBox , WindFarm):
    """
    Note: this sizing imply a choice of the reference wich is obvious if all turbines are same => there is ONE Radius
    :param MannBox:
    :param R_WindTurbine:
    :param U_mean_WindFarm:
    :param windFarm_lenght:
    :return:
    """

    MannBox.U = MannBox.U / WindFarm.U_mean
    print 'U (no dimension): ', MannBox.U
    MannBox.L = MannBox.L / WindFarm.WT_R
    print 'L (no dimension): ', MannBox.L

    MannBox.T = MannBox.L / MannBox.U; print 'Mannbox T_total (no dimension): ', MannBox.T
    MannBox.dt = MannBox.T / MannBox.nx; print 'dt (no dimension): ', MannBox.dt
    MannBox.ti = [i*MannBox.dt for i in range(MannBox.nx)]
    MannBox.ti = [t for t in MannBox.ti
                  if t < MannBox.SimulationTime + MannBox.dt]  # we catch just one point after the simulation time


    # Should depends on the FFoR domain.
    MannBox.lx = MannBox.L
    MannBox.ly = MannBox.ly / WindFarm.WT_R * k
    MannBox.lz = MannBox.lz / WindFarm.WT_R * k

    MannBox.dx = MannBox.lx / MannBox.nx
    MannBox.dt = MannBox.dx / MannBox.U

    MannBox.u_TurbBox = MannBox.u_TurbBox / MannBox.U_ref
    return MannBox






