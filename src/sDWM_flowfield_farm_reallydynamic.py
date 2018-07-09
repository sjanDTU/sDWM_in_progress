"""
author: Augustin Grosdidier
date: 06/06/2018
Purpose:
from the original sDWM code (steady-state), develop a real Dynamic sDWM taking account of the meandering with Mann box
and wake added turbulence.

Modified from the Main DWM core program: flowfield calculation
created by Ewan Machefaux <ewan.machefaux@gmail.com>
"""

import numpy as np
import matplotlib.pyplot as plt

def DWM_make_grid(meta):
    """
    Originally present in DWM_flowfield_farm developped by Ewan Machefaux <ewan.machefaux@gmail.com>:

    Function to adapt the grid (polar MFoR and cartesian FFoR) length
    in the streamwise direction based on the distance to the next turbine
    It also create extraction plane in the cartesian ffor domain

    Updated:
    - take account dynamic physic, therefore we add a time-dimension.
    Parameters
    ----------
    meta    class holding grid and ambient parameters

    Returns
    -------
    meta (updated)
            lz_mixl  : length of the mixing length domain
            vz_mixl  : vector of streamwise grid points of the polar MFoR domain in R
            vz       : vector of streamwise grid points indices of the cartesian domain
            nz       : vector length of streamwise grid points of the cartesian domain in D
            z_vec    : vector length of streamwise grid points of the cartesian domain in D
    """

    meta.vz=meta.hub_z[0:]*meta.dz
    #print 'meta.vz='+str(meta.vz)
    # The mixL domain
    meta.lz_mixl=1.0+2.0*(max(meta.hub_z)) # in R: 1R longer than ffor flow field due to backward diff scheme
    meta.vz_mixl = np.linspace(0,meta.lz_mixl-meta.dz_mixl,meta.lz_mixl/meta.dz_mixl) # coordinate streamwise in R vector mixL
    ### Plotting
    #print "dz_mixl="+str(meta.dz_mixl)
    #print "vz_mixl="+str(meta.vz_mixl)
    #print "len(vz_mixl)="+str(len(meta.vz_mixl))

    meta.nz = len(meta.vz)           # nb points in z (streamwise) direction ffor flow field
    meta.z_vec= meta.vz/meta.dz

    return meta

def DWM_MFOR_to_FFOR(mfor,meta,meand,ffor):
    """
    Originally present in DWM_flowfield_farm developped by Ewan Machefaux <ewan.machefaux@gmail.com>:
    Function that calculate the velocity in the fixed (global) frame of reference from the Mfor

    Updating:
    taking account of the meandering effect due to a mann box
    Parameters
    ----------
    meand (instance of class)   class holding the meandering parameters
    meta (instance of class)    class holding grid and ambient parameters

    Returns
    -------
    ffor (instance of class)    updated class holding global velocity field
        x_vec(1,nx) : coordinates in X direction
        y_vec(1,ny) : coordinates in Y direction
        z_vec(1,nz) : coordinates in Z direction (streamwise)
        x_mat(nx,ny): coordinates matrix in meshgrid format for the X component
        y_mat(nx,ny): coordinates matrix in meshgrid format for the Y component
        z_mat(nx,nz): coordinates matrix in meshgrid format for the Z component
        TI_meand_axial_ffor (nx,ny,nz): turbulence due to wake meandering in global coordinate system
        WS_axial_ffor (nx,ny,nz): velocity deficit in global coordinate system
        TI_axial_ffor (nx,ny,nz): apparent turbulence in global coordinate system see Madsen et al [2]
    meta (instance of class)
    """

    ##############################################################################################################
    # recalculate into Cartesian grid
    # initiate/reset Cartesian flow field
    print 'Performing MFoR to FFoR Computation'
    ffor.ffor_flow_field_TI_tmp_tmp =  meta.TI * np.ones((meta.nx,meta.ny))  #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.TI_axial_ffor_tmp     =  np.zeros((meta.nx,meta.ny,meta.nz))        #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.WS_axial_ffor_tmp     =  np.zeros((meta.nx,meta.ny,meta.nz))        #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.ffor_flow_field_ws_tmp2    =  np.zeros((meta.nx,meta.ny,meta.nz))   #X = lateral ,Y = vertical, time ,Z = streamwise
    print 'shape(ffor.WS_axial_ffor_tmp): ', np.shape(ffor.WS_axial_ffor_tmp)
    ffor.TI_meand_axial_ffor    =  np.zeros((meta.nx,meta.ny,meta.nz))
    ffor.WS_axial_ffor =  np.zeros((meta.nx,meta.ny,meta.nz))
    ffor.TI_axial_ffor =  np.zeros((meta.nx,meta.ny,meta.nz))
    print 'shape(ffor.WS_axial_ffor): ', np.shape(ffor.WS_axial_ffor)
    ffor.x_vec_t            =np.zeros((meta.nx,meta.nz))
    ffor.x_mat_t            =np.zeros((meta.nx,meta.ny,meta.nz))
    # Define radial distance vectors
    # r_dist_2              = 999*np.ones((meta.nx,meta.ny))
    r_dist                = 999*np.ones((meta.nx,meta.ny))

    # CREATES THE GLOBAL FLOW FIELDS IN CARTESIAN GRID
    # meta.x_mat = np.tile(meta.x_vec.reshape(len(meta.x_vec),1),meta.ny)
    # meta.y_mat = np.tile(meta.y_vec,(meta.nx,1))
    meta.z_mat = np.tile(meta.vz,(meta.nx,1))

    for i_z in np.arange(0,meta.nz,1):
        # EXTRACT TI_DWM AND WS_DWM IN MFoR
        #Plot deficit MFOR
        try:
            DWM_WS_DATA = mfor.U[meta.vz[i_z],:]

            print 'DWM_WS_DATA: ', DWM_WS_DATA


        except:
            print 'Fatal error, possibly due to a too low turbulence intensity with respect to the demanded mean wind speed, try increase the input TI'
        DWM_TI_DATA = mfor.TI_DWM[meta.vz[i_z],:]
        print 'DWM_TI_DATA: ', DWM_TI_DATA

        #Plot TI MFOR

        """
        plt.plot(meta.vr_mixl,mfor.TI_DWM[meta.vz[i_z], :])
        plt.title('TI in MFoR at Turbine '+str(7-i_z)+' Location')
        plt.xlabel('vr (polar discretization)')
        plt.ylabel('TI (Turbulence Intensity in MFoR)')
        plt.show()
        """

        ### Correct DWM_TI_DATA so that no point to have lower TI than "TIamb"
        DWM_TI_DATA[DWM_TI_DATA < np.nanmean(meta.mean_TI_DWM)] = np.nanmean(meta.mean_TI_DWM)
        plt.ion()
        plt.figure()
        plt.title('Center?')

        # Change this loop to be more suitable with the pseudo-lagragian approach described in the pragmatic approach
        # of the wake meandering. It involve to change the meand data wich is originally based on DWM_meta_meand
        for i_t in np.arange(0,len(meand.time),1):

            Ro_x                    = meand.meand_pos_x[i_z,i_t]
            Ro_y                    = meand.meand_pos_y[i_z,i_t]

            print '(Ro_x,Ro_y): ', [Ro_x, Ro_y]
            plt.plot(Ro_x,Ro_y,'x')
            plt.draw
            plt.pause(0.1)

            r_dist                  = np.sqrt((meta.x_mat - Ro_x)**2 + (meta.y_mat - Ro_y)**2 )
            print 'r_dist: ', r_dist

            tmp_index               = r_dist < mfor.WakeW[meta.vz[i_z]]*1.5
            tmp_field_WS            = np.ones((meta.nx,meta.ny))
            print 'tmp_index: ', tmp_index

            # it's here that we change the velocity to be in FFOR
            tmp_field_WS[tmp_index] = np.interp( r_dist[tmp_index],meta.vr_m, DWM_WS_DATA)
            ffor.WS_axial_ffor_tmp[:, :, i_z]  = ffor.WS_axial_ffor_tmp[:, :, i_z] + (tmp_field_WS)
            ffor.ffor_flow_field_ws_tmp2[:, :, i_z] = ffor.ffor_flow_field_ws_tmp2[:, :,i_z]+ (tmp_field_WS**2)

            tmp_field_TI            = meta.TI * np.ones((meta.nx,meta.ny))
            tmp_field_TI[tmp_index] = np.interp( r_dist[tmp_index],meta.vr_m,DWM_TI_DATA)

            ffor.ffor_flow_field_TI_tmp_tmp[:, :]      = tmp_field_TI
            ffor.TI_axial_ffor_tmp[:, :, i_z]     = ffor.TI_axial_ffor_tmp[:, :, i_z] + ffor.ffor_flow_field_TI_tmp_tmp**2
        plt.ioff()
    # Stores the mean field
    for i_z in np.arange(0,meta.nz,1):
        #### Here we average all time related data along the z-axis
        #### That not for a Dynamic Expression, we want to store the field depending on time to be a good input for
        #### flex5 Computation
        ffor.TI_meand_axial_ffor[:, :, i_z]=np.sqrt(abs(ffor.ffor_flow_field_ws_tmp2[:, :, i_z] - ((ffor.WS_axial_ffor_tmp[:, :, i_z]**2)/len(meand.time)) )/ (len(meand.time)-1.0))
        ffor.WS_axial_ffor[:, :, i_z]      = (ffor.WS_axial_ffor_tmp[:, :, i_z]  / len(meand.time))
        ffor.TI_axial_ffor[:, :, i_z]      = (ffor.TI_axial_ffor_tmp[:, :, i_z]  / len(meand.time))**(1.0/2.0)
        ffor.x_vec_t[:, i_z]               = (meta.x_vec-meta.hub_x[i_z])/2.
        ffor.x_mat_t[:,:,i_z]              = np.tile(ffor.x_vec_t[:, i_z] .reshape(len(ffor.x_vec_t[:, i_z]),1),meta.ny)/2.

    # Store the ffor flow field
    ffor.x_vec                   = (meta.x_vec-meta.hub_x[0])/2.
    ffor.y_vec                   = (meta.y_vec-meta.hub_y)/2.
    ffor.z_vec                   = meta.z_vec+np.hstack((0., np.cumsum(meta.hub_z[0:])))[0]
    ffor.x_mat                   = meta.x_mat/2.
    ffor.y_mat                   = meta.y_mat/2.
    ffor.z_mat                   = meta.z_mat/meta.dz

    return mfor,ffor,meta,meand