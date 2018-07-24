"""
author: Augustin Grosdidier
date: 02/07/2018

Running Core of the meandering computation process.
Computing for MannBox and at the end for LESBox.

Mainly Based on 'A Pragmatic Approach oh the Meandering' GLarsen
Work with the 'Simplification' for the moment (very fast and relevant for a MannBox)
But it could be improved with no simplification easily.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import math as m
from matplotlib import animation

from DWM_GClarsenPicks import get_Rw

import WindTurbine as wt

from ReadTurbulence import *

from cTurb import *

from Polar_Interpolate_Data_Plan import *

import time

def DWM_extract_meandering_from_TurbBox(MannBox, WindFarm):
    """
    Notice: in Future Work we can add the fact that the Turbulence Box becomes slower through the windfarm
    :param filename:
    :param WindFarm:
    :param WT:
    :return:
    """
    # WT location in the streamwise (come from meta.z_vec of Make Grid in sDWM):
    # [0  4  8 12 17 21 25 30]
    # ---- # For Test, this have to come from sDWM and so to be a parameter # ---- #
    WindFarm.stream_location_z = [2*ld for ld in WindFarm.stream_location_z] # [R]
    # reverse index location: originally the wind farm is index like this: 7 6 5 4 3 2 1 0 (upstream to downstream)
    WindFarm.stream_location_z = WindFarm.stream_location_z[::-1]
    WindFarm.stream_location_z = [abs(ld-max(WindFarm.stream_location_z)) for ld in WindFarm.stream_location_z]
    #                       : now 0 1 2 3 4 5 6 7
    print 'WindFarm.stream_location_z ', WindFarm.stream_location_z
    #then [0,  4,  8, 12, 17, 21, 25]
    # then [0,  4,  8, 12, 17, 21]
    # then [0,  4,  8, 12, 17]
    # then [0,  4,  8, 12] etc...
    WindFarm.nodim_lenght = WindFarm.stream_location_z[-1]

    video = False

    # -------------------------------- # INITIALIZATION # ------------------------------------ #
    MannBox, WindFarm = Init_Turb(MannBox, WindFarm)

    # -------------------------- # PLOT CHECK PART # ---------------------------------- #
    # With Simplification described in 'Pragmatic approach of wake meandering'
    """
    Ld = 100./WindFarm.WT_R
    ti = 1 #s
    yc_zc = Wake_Dynamic_at_Ld(Ld, ti, MannBox, Meand_Mann)
    print yc_zc
    yc_zc = []
    for ti in MannBox.ti[:200]:
        print 'ti: ', ti
        yc_zc.append(Wake_Dynamic_at_Ld(Ld, ti, MannBox, Meand_Mann))
    yc_zc = np.array(yc_zc)
    print yc_zc
    print np.shape(yc_zc)  # first axis: time, second axis: [vc, wc]
    plt.figure()
    plt.plot(yc_zc[:, 0],yc_zc[:, 1], label='Center location')
    plt.show()
    #"""

    # do a check of the 'Perfect Correlation' described in 'Pragmatic Approach

    #### Plot yc or zc in time in WindFarm (3D) with circle of 2*Rw radius
    """
    Ld = np.linspace(0, WindFarm.lenght/WindFarm.WT_R, 200)
    Wake = []
    T = [ ti for ti in MannBox.ti[:] if ti<MannBox.SimulationTime]

    for ld in Ld:
        print 'Ld: ', ld
        yc_zc = []
        # as x_b = U*tho, tho = x_b/U (here x_b = ld)
        tho = ld / MannBox.U    # time for the first wake release to reach the studied plan at ld

        # Get Wake radius function of the transportation time (seems at this distance)
        if MannBox.WakeExpansion:
            # ----# Get the wake radius #----#
            Rw = get_Rw(x=ld, R=1., TI=MannBox.TI, CT=WindFarm.CT)  # resulted Rw is no dimensionalized
            MannBox.WakeRadius_for_this_Plan = Rw
            MannBox.Af = m.pi * (2 * MannBox.WakeRadius_for_this_Plan) ** 2

        for ti in T:
            #print 'ti: ', ti
            #print 'tho: ', tho

            # a way to do well is to make difference between ti and tho
            #ti = ti - tho
            #print 'ti: ', ti
            if ti-tho < 0:
                yc_zc.append([np.nan, np.nan])
            else :
                yc_zc.append(Wake_Dynamic_at_Ld(ld, ti-tho, MannBox, Meand_Mann, WindFarm))
            # a way to do well is to make difference between ti and tho
            ti = ti - tho
        Wake.append(yc_zc)
    Wake = np.array(Wake)
    print 'shape Wake', np.shape(Wake)

    # Compute wake radius for each distances
    L_Total = np.array([MannBox.dx*i for i in range(0, MannBox.nx, 64) if i*MannBox.dx < WindFarm.lenght/WindFarm.WT_R])
    twoRw_to_plot = 2*get_Rw(L_Total, R=1., TI=MannBox.TI, CT=WindFarm.CT)
    y_rw, z_rw = [], []
    for index_ld in range(len([MannBox.dx*i for i in range(0, MannBox.nx, 64) if i*MannBox.dx < WindFarm.lenght/WindFarm.WT_R])):
        y_rw.append(twoRw_to_plot[index_ld] * np.cos(np.linspace(-m.pi, m.pi,15)))
        z_rw.append(twoRw_to_plot[index_ld] * np.sin(np.linspace(-m.pi, m.pi, 15)))
    y_rw, z_rw = np.array(y_rw), np.array(z_rw)
    print 'shape y_rw: ', np.shape(y_rw)

    raw_input('Press any key to continue and plot solutions: ')
    plt.ion()
    fig = plt.figure()

    for index_t in range(len(T)):
        if video:
            break
        yc_zc = Wake[:, index_t, :]
        print 'index_t', index_t
        if index_t == 2 or index_t == 50 or index_t == 200 or index_t == 210:
            plt.pause(10)

        plt.pause(0.1)
        plt.clf()
        plt.cla()
        ax = fig.add_subplot(111, projection='3d')
        view = MannBox.ly/2
        plt.ylim(-view, view)
        plt.xlim(-view, view)

        ax.scatter(yc_zc[:, 0], yc_zc[:, 1], Ld)
        #ax.scatter for circle of 2 wake radius centered with the wake center
        for i in range(np.shape(y_rw)[1]):
            ax.scatter(y_rw[:, i] , z_rw[:, i], L_Total, color='r')
        ax.view_init(azim=0, elev=-90)
        #ax.set_zlim(0, WindFarm.lenght/WindFarm.WT_R/10)
        ax.set_zlim(0, WindFarm.lenght/WindFarm.WT_R)
        plt.draw()
        #plt.pause(1)
    #"""

    #### Plot vc or wc in time at a specified distance (account wake expansion)
    """
    ld = 4   #(4R)
    #ld = 100 / WindFarm.WT_R  # We are at a specified distance so we can get the wake radius at this distance
                              # (doesn't change in time)
                              # A simple stationary semi-analytical wake model Gunner C. Larsen (2009)

    T = [ti for ti in MannBox.ti[:] if ti < MannBox.SimulationTime]

    print 'Ld: ', ld
    vc_wc = []
    # as x_b = U*tho, tho = x_b/U (here x_b = ld)
    tho = ld / MannBox.U  # time for the first wake release to reach the studied plan at ld

    # Get Wake radius function of the transportation time (seems at this distance)
    if MannBox.WakeExpansion:
        # ----# Get the wake radius #----#
        # resulted Rw is no dimensionalized
        MannBox.WakeRadius_for_this_Plan = get_Rw(x=ld, R=1., TI=MannBox.TI, CT=WindFarm.CT)
        MannBox.Af = m.pi * (2 * MannBox.WakeRadius_for_this_Plan) ** 2

    for ti in T:
        print 'ti: ', ti
        if ti - tho < 0:
            vc_wc.append([ti, np.nan, np.nan])
        else:
            vc_wc.append([ti]+(Wake_Dynamic_at_Ld(ld, ti - tho, MannBox, Meand_Mann, WindFarm)))

    vc_wc = np.array(vc_wc)
    yc_zc = vc_wc
    yc_zc[:, 1:3] = ld/MannBox.U * yc_zc[:, 1:3]

    print 'vc_wc shape: ', np.shape(vc_wc)

    print 'yc_zc: '
    print yc_zc

    plt.figure()
    plt.title('Characteristic Velocities in Time at a specified distance')
    plt.plot(vc_wc[:, 0], vc_wc[:, 1], label='vc(t) at ld='+str(ld))
    plt.plot(vc_wc[:, 0], vc_wc[:, 2], label='wc(t) at ld='+str(ld))
    plt.legend()
    plt.show()

    plt.figure()
    plt.title('center movements in time at a specified distance')
    plt.plot(yc_zc[:, 0], yc_zc[:, 1], label='yc(t) at ld='+str(ld))
    plt.plot(yc_zc[:, 0], yc_zc[:, 2], label='zc(t) at ld=' + str(ld))
    plt.plot(yc_zc[:, 0], np.sqrt(yc_zc[:, 1]**2 + yc_zc[:, 2]**2), label='rc(t) at ld=' + str(ld))
    plt.legend()
    plt.show()
    #"""

    # -------------- # Meandering Computation for each plan of interest # --------------------------- #
    # Ld is the list of the distance between generating plan and the plans of interest
    #Ld = WindFarm.stream_location_z  # We are at a specified distance so we can get the wake radius at this distance
                                     # (doesn't change in time)
                                     # A simple stationary semi-analytical wake model Gunner C. Larsen (2009)
    Meand_Mann = meand_mann()
    if MannBox.CorrectionDelay:
        # We want to begin the simulation when the first plan go out of the WindFarm Box
        delay = WindFarm.nodim_lenght / MannBox.U
    # The data do not exceed the simulation time defined in class object (cMeand or cMann)
    T = [ti for ti in MannBox.ti[:] if ti < MannBox.SimulationTime]


    WAKES = [] # list holding all wakes generated by each turbines
    for WT_index in range(len(WindFarm.stream_location_z)):
        # creates the list ld of interest
        # all the distances are calculate for reference WT, we delete WT upstream to the WT ref
        ld_ref = WindFarm.stream_location_z[WT_index]
        Ld = [abs(ld-ld_ref) for ld in WindFarm.stream_location_z[WT_index:]]

        # For each ld in Ld, we get a data matrix structured like this: column 1: time(s), column 2: yc, column 3: zc
        # we store these matrixs in a list along Ld

        Wake = []

        for ld in Ld:
            print 'Ld: ', ld

            ts_vc_wc =[]
            # as x_b = U*tho, tho = x_b/U (here x_b = ld)
            tho = ld / MannBox.U  # time for the first wake release to reach the studied plan at ld

            # Get Wake radius function of the transportation time (seems at this distance)
            if MannBox.WakeExpansion:
                # ----# Get the wake radius #----#
                Rw = get_Rw(x=ld, R=1., TI=MannBox.TI, CT=WindFarm.CT)  # resulted Rw is no dimensionalized
                MannBox.WakeRadius_for_this_Plan = Rw
                MannBox.Af = m.pi * (2 * MannBox.WakeRadius_for_this_Plan) ** 2
            boolref = False
            for ts in T:
                if MannBox.CorrectionDelay:
                    ts = ts + delay
                    vc_wc = Wake_Dynamic_at_Ld(ld, ts - tho, MannBox, Meand_Mann)
                    ts_vc_wc.append([ts-delay] + vc_wc)
                else:
                    if ts - tho < 0:
                        ts_vc_wc.append([ts, np.nan, np.nan])
                    else:

                        vc_wc = Wake_Dynamic_at_Ld(ld, ts - tho, MannBox, Meand_Mann, WindFarm)
                        if not boolref:
                            boolref = True
                            vc_wc_ref = vc_wc
                        ts_vc_wc.append([ts]+vc_wc)

            ts_vc_wc = np.array(ts_vc_wc)
            ts_yc_zc = ts_vc_wc
            if not MannBox.CorrectionDelay:
                NaN_index = np.isnan(ts_yc_zc[:, 1])
                ts_yc_zc[NaN_index, 1] = vc_wc_ref[0]
                ts_yc_zc[NaN_index, 2] = vc_wc_ref[1]

            ts_yc_zc[:, 1:3] = ld/MannBox.U * ts_yc_zc[:, 1:3]
            Wake.append(ts_yc_zc)
        WAKES.append(Wake)

    # Index of the list WAKES correspond directly of the iteration in sDWM

    #Plot Part
    # Plot each Wake
    if MannBox.RESULT_plot:
        Plot_each_Wake(WAKES, WindFarm)
        Plot_Wakes_at_each_ld(WAKES, WindFarm)
    DATA = WAKES
    np.save('WAKES',DATA)
    print 'Wake Radius Data saved...'
    return

def Init_Turb(MannBox, WindFarm):

    video = False

    # -------------------------- # SIZE THE TURBULENT BOXES # ------------------------ #
    # sizing mann for the biggest wake?
    # Get the max wake radius to size correctly the Mann box,
    # (be careful if you don't use the 'simplification' described by Larsen,
    # you must manage the wake movement in addition of the wake radius to size the turbulent box.
    # (note: not yet implemented but we can also positionate rightly the turbulent Box to be close to the ground
    # and we could manage the ground effect on the meandering)
    # for the moment the wake movement doesn't need to take care of the ground (no need to implement reflecting surface)
    WindFarm.Rw_Max = get_Rw(x=WindFarm.nodim_lenght, R=1., TI=MannBox.TI, CT=WindFarm.CT)
    print 'Rw max in WindFarm is: ', WindFarm.Rw_Max

    MannBox = sizing_MannBox(MannBox, WindFarm)  # put: , windFarm_lenght=0.) for the first function of  Main loop
    print 'ly/2: ', MannBox.ly/2
    print 'dt: ', MannBox.dt

    # -------------------------- # INITIALIZATION # --------------------------------- #

    print 'Total simulation time for the entire MannBox: ', MannBox.ti[-1]  # 819,2s in total

    # ---- # Get the wake radius # ---- #
    MannBox.WakeRadius_for_this_Plan = 1. # the wake generating plan rotor radius = wake radius (we assumed)
    MannBox.Af = m.pi * (2 * MannBox.WakeRadius_for_this_Plan) ** 2

    ################################################################################################

    # Initialization if no 'Simplification'
    """
    # ---- # Compute for the two component for the generating plan # ---- #
    # According to 'Wake Meandering: A Pragmatic Approach' (2008) by Larsen and co.
    # We deal just with v, w component
    for char in ['vfluct', 'wfluct']:
        MannBox.plan_of_interest = get_plan_of_interest(MannBox, ti=0, component_char=char)

        # ---- # Interpolation Part (WakeCentered) # ---- #
        Interpo_Integrate = interpo_integrate()
        Interpo_Integrate = Interpolate_plan(MannBox, Interpo_Integrate, Meand_Mann)

        # ---- # Integration Part # ---- #
        Interpo_Integrate = polar_mesh_value(MannBox, Interpo_Integrate)
        Final_Integral_Value = Trapz_for_Integrate_general_grid(Interpo_Integrate)

        # ---- # Calculate the characteristic velocity # ---- #
        Meand_Mann.init_vc_wc.append(Final_Integral_Value / MannBox.Af)  # [vc, wc]
    #"""

    return MannBox, WindFarm

def Wake_Dynamic_at_Ld(ld, ti, MannBox, Meand_Mann):
    """
    Based on Simplification part of " Pragmatic Approach of wake meandering "
    Double Assumptions:
        - Restricted to wake deficit behaviour at a given downstream distance from the wake generating Wind Turbine
        - Perfect correlation between the characteristic transversal and vertical velocities in all 'cross section'
        of the relevant turbulent box /!\

    "
    A priori, the simplification is believed to be a reasonable approximation for moderate downstream distances,
    where the wake centre displacements are modest, and where the change in the driving large-scale turbulence
    structures therefore also is moderate in the relevant spatial regime. For larger downstream distances, account
    should be taken to the spatial variability of the large-scale turbulence components, and the simplification
    consequently breaks down.
    "
    :param Ld: float
    :return:
    """
    # Caluculation to get Center position directly
    bool_center = False
    # Calculation to get Characteristic Velocities (this way have to be used for the MEANDERING MAIN)
    bool_velocity = True
    # (Ld, y_g(Ld/U; t0+tho), z_g(Ld/U; t0+tho)
    # (y_g, z_g) = Ld/U * [vc(U(T-ti), 0, 0), wc(U(T-ti),0,0)]


    # Determine vc and wc at U(T-ti), 0, 0:
    yc_zc = []
    vc_wc = []
    for char in ['vfluct', 'wfluct']:
        Meand_Mann.wake_center_location = (0, 0) # due to simplification
        MannBox.plan_of_interest = get_plan_of_interest(MannBox, ti=ti, component_char=char)

        # ----# Interpolation Part (WakeCentered) #----#
        Interpo_Integrate = interpo_integrate()
        Interpo_Integrate = Interpolate_plan(MannBox, Interpo_Integrate, Meand_Mann)

        # ----# Integration Part #----#
        Interpo_Integrate = polar_mesh_value(MannBox, Interpo_Integrate)
        Final_Integral_Value = Trapz_for_Integrate_general_grid(Interpo_Integrate)

        # ----# Calculate the characteristic velocity #----#
        if bool_center:
            yc_zc.append(ld/MannBox.U * Final_Integral_Value / MannBox.Af)
        if bool_velocity:
            vc_wc.append(Final_Integral_Value / MannBox.Af)
    if bool_center:
        return yc_zc
    if bool_velocity:
        return vc_wc
########################################################################################################################

#"""
