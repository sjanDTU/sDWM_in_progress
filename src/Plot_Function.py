"""
author: Augustin Grosdidier
date: 11/07/2018
some usual plot functions for Meandering
"""
import numpy as np
import matplotlib.pyplot as plt
import math as m
def plot_interpolation(MannBox, Interpo_Integrate):
    colorbar_boundaries_setting = 0.3
    cbs = colorbar_boundaries_setting
    val_mat_pol = Interpo_Integrate.f_cart(Interpo_Integrate.CartMesh[0], Interpo_Integrate.CartMesh[1])
    ypp, zpp = Interpo_Integrate.CartMesh

    fig = plt.figure()
    ax = fig.add_subplot(122)
    plt.contourf(ypp, zpp, val_mat_pol, np.arange(-cbs, cbs, .01), cmap='jet')
    plt.colorbar()
    ax.set_xlabel('y'), ax.set_ylabel('z')

    y = np.linspace(-MannBox.ly / 2, MannBox.ly / 2, MannBox.ny)
    z = np.linspace(-MannBox.lz / 2, MannBox.lz / 2, MannBox.nz)
    yy, zz = np.meshgrid(y, z)

    plt.subplot(121)
    plt.contourf(yy, zz, MannBox.plan_of_interest, np.arange(-cbs, cbs, .01), cmap='jet')
    plt.plot(2*MannBox.WakeRadius_for_this_Plan * np.cos(np.linspace(-m.pi, m.pi)),
             2*MannBox.WakeRadius_for_this_Plan * np.sin(np.linspace(-m.pi, m.pi)), 'k', label='Af aera')
    plt.colorbar()
    ax.set_xlabel('y'), ax.set_ylabel('z')

    plt.show()
    return

def Plot_each_Wake(WAKES, WindFarm):
    for i_wake in range(len(WAKES)):
        plt.figure(i_wake)
        plt.title('Lateral position for the Wake ' + str(i_wake))
        plt.xlabel('Time [s]')
        plt.ylabel('y [R]')
        for i_ld in range(len(WAKES[i_wake])):
            plt.plot(WAKES[i_wake][i_ld][:, 0], WAKES[i_wake][i_ld][:, 1], label='at WT'+str(i_ld+i_wake))
        plt.legend()
    plt.show()
    return

def Plot_Wakes_at_each_ld(WAKES, WindFarm):
    for i_ld in range(len(WindFarm.stream_location_z)):
        print 'i_ld:', i_ld
        plt.figure(i_ld)
        plt.title('Lateral position for each Wake at ld = ' + str(WindFarm.stream_location_z[i_ld]))
        plt.xlabel('Time [s]')
        plt.ylabel('y [R]')
        for i_wake in range(i_ld+1):
            print 'i_wake ', i_wake
            print 'i_ld ', i_ld
            #i_wake = i_wake - i_ld
            plt.plot(WAKES[i_wake][i_ld][:, 0], WAKES[i_wake][i_ld][:, 1],label='Wake ' + str(i_wake))
            i_ld = i_ld - 1
        plt.legend()
    plt.show()