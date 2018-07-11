# -*- coding: utf-8 -*-
""" Main DWM core program: flowfield calculation
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""

import numpy as np
import time
from cDWM import Meta, Aero, Meand, FFoR, MFoR, Outputs
from cBEM import InitBEM
from DWM_calc_mixL import DWM_calc_mixL

# Choose your Bem calculation code
#from DWM_main_BEM import getInduction
from DWM_main_BEM_new_controller import getInduction


from math import pi, sqrt, isnan
from scipy import io, interpolate
from scipy.integrate import simps
from DWM_misc import smooth

import matplotlib.pyplot as plt


def DWM_main_field_model(ID_waked,deficits,inlets_ffor,inlets_ffor_deficits,inlets_ffor_turb,turb,DWM,out,**par):
    """Main flow field calculation function, handling all the calls to each sub functions. This function is called for
       each turbine in the wind farm from the most upstream to the most downstream one. The flow field calculations are
       carried out at the downstream distance of interest, i.e., where downstream rotors are partially or fully in the
       wake of the upstream turbine.

        Inputs
        ----------
        ID_waked: dict(nWT) holding list of upstream turbine index for each turbine in the wind farm
        deficits: dict(nWT) holding a list of the mean rotor averaged deficit from upstream wakes
        turb:  dict(nWT) holding a list of mean turbulence intensities contributions from upstream wakes
        inlets_ffor: dict(nWT) holding a list of array containing the flow field in the fixed frame of reference from
        upstream wakes contributions
        inlets_ffor_deficits: dict(nWT) holding a list of array containing the flow field in the fixed frame of
        reference from upstream wakes contributions at the rotor position
        inlets_ffor_turb: dict(nWT) holding a list of array containing the turbulence field in the fixed frame of
        reference from upstream wakes contributions at the rotor position
        out: dict(nWT),holding the main sDWM model outputs i.e mean power from BEM, mean power estimated from powercurve,
        mean rotor averaged wind speed, mean rotor average turbulence intensity, mean thrust coefficient from BEM and
        from power curve
        par: dict(nWT) holding the DWM cylindrical and cartesian grid coordinates as well as turbine and ambient conditions

        Outputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        aero (instance of class): Instance of class Aero holding BEM-aero core variables
        mfor (instance of class): Instance of class Mfor holding the meandering frame of reference scalars used by the
        Ainslie model in local Mixl coordinates
        ffor:(instance of class): Instance of class Ffor holding the fixed frame of reference velocity field in global
        WF coordinates
        DWM: dict(nWT) list containing full outputs of the sDWM (including flow field in ffor and mfor) See description
        of DWM_outputs for more details
        deficits: dict(nWT) update list of deficits contributions from upstream wakes
        inlets_ffor: dict(nWT) updated list of array containing the flow field in the fixed frame of reference from
        upstream wakes contributions
        inlets_ffor_deficits: dict(nWT) updated list of array containing the flow field in the fixed frame of reference
        from upstream wakes contributions at the rotor position
        inlets_ffor_turb: dict(nWT) updated list of array containing the turbulence field in the fixed frame of
        reference from upstream wakes contributions at the rotor position
        turb:dict(nWT) updated list of mean turbulence intensities contributions from upstream wakes
        out: returns the mean power from BEM, mean power estimated from power curve,mean rotor averaged wind speed, mean
        rotor average turbulence intensity,
        mean thrust coefficient from BEM and from power curve
        ID_waked: dict(nWT) holding list of upstream turbine index for each turbine in the wind farm
    """
    # Create instances of class
    meta=Meta()
    meta.parse(**par)
    meand=Meand()
    ffor =  FFoR()
    aero = Aero(meta.WTG)
    t = time.time()
    # ##### Set up MFoR and FFoR streamwise domain properties   #####################################################
    meta                 = DWM_make_grid(meta)
    # ##### Load wake meandering properties from meta model: f(stab,hub height,z,TI) ################################
    if meta.previous_sDWM:
        meand                = DWM_meta_meand(meand,meta)
    else:
        meta, meand = get_Meandering_dynamic(meta, meand)

    # ##  Run BEM model and create velocity deficit calculations inputs #############################################
    start_time = time.time()
    if meta.previous_sDWM:
        aero,mfor,out,BEM    = DWM_aero(meta,ffor,aero,deficits,turb,inlets_ffor,inlets_ffor_deficits,out,ID_waked)
    else:
        if meta.steadyBEM_AINSLIE:
            aero, mfor, out, BEM = DWM_aero(meta, ffor, aero, deficits, turb, inlets_ffor,inlets_ffor_deficits,out,ID_waked)
        elif not meta.steadyBEM_AINSLIE:
            aero, mfor, out, BEM = DWM_aero_dynamic(meta, ffor, aero, deficits, turb, inlets_ffor, inlets_ffor_deficits, out, ID_waked)
    print 'Computation Time for BEM is: ', time.time() - start_time

    # ############## Perform wake velocity calculations in MFoR #####################################################

    start_time = time.time()
    mfor                 = DWM_calc_mixL(meta,aero,mfor)
    print 'Computation Time for Ainslie is: ', time.time()-start_time
    # ############# Reconstruct global flow field by applying wake meandering #######################################
    if meta.previous_sDWM:
        mfor,ffor,meta,meand = DWM_MFOR_to_FFOR(mfor,meta,meand,ffor)
    else:
        mfor, ffor, meta, meand = DWM_MFOR_to_FFOR_dynamic(mfor, meta, meand, ffor)

    # ############## Compute deficit at downstream rotor
    if meta.previous_sDWM:
        deficits, ID_waked,inlets_ffor,inlets_ffor_deficits   = DWM_get_deficit(ffor,meta,deficits,ID_waked,inlets_ffor,inlets_ffor_deficits)
    else:
        deficits, ID_waked, inlets_ffor, inlets_ffor_deficits  = DWM_get_deficit_FFOR_dynamic(ffor, meta, deficits, ID_waked, inlets_ffor, inlets_ffor_deficits)

    # print deficits
    # print inlets_ffor
    # ############## Compute turbulence level at downstream rotor

    if meta.previous_sDWM:
        turb, inlets_ffor_turb                 = DWM_get_turb(ffor,meta,turb,inlets_ffor_turb)
    else:
        turb, inlets_ffor_turb = DWM_get_turb_dynamic(ffor, meta, turb, inlets_ffor_turb)

    # ############# Compute new axisymmetric flow for next turbine ##################################################
    # ffor,inlets = DWM_make_inflow_to_mixl(meta,ffor,inlets)
    # ############# Write relevant results in DWM variables #########################################################
    # if meta.full_output is True:
    # DWM                  = DWM_outputs(DWM,ffor,mfor,meta,aero,par,BEM)
    elapsed = time.time() - t
    print '*********Turbine %i (%i turbine in its wake) produces %4.2f kW at %4.2f m/s completed in %4.2f sec ***********************************' %(meta.wtg_ind[0],len(meta.wtg_ind[1:]),aero.pow_cur,meta.mean_WS_DWM,elapsed)
    print '****Turbine %i (%i turbine in its wake) produces %4.2f kW at %4.2f m/s completed in %4.2f sec **********' %(meta.wtg_ind[0],len(meta.wtg_ind[1:]),aero.Power/1000.,meta.mean_WS_DWM,elapsed)
    return( aero, meta, mfor, ffor, DWM, deficits,inlets_ffor, inlets_ffor_deficits,inlets_ffor_turb,turb, out,ID_waked)


########################################################################################################################
# ***********************************************Previous sDWM**********************************************************
# *******************************Working with statistical approach of meandering****************************************
########################################################################################################################

def DWM_aero(meta,ffor,aero,deficits,turb,inlets_ffor,inlets_ffor_deficits,out,ID_waked):
    """ Aerodynamique module of the DWM. This module contains the wake summation module (deficit and turbulence accumulation)
    The steady state blade element momentum

        Inputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        aero (instance of class): Instance of class Aero holding BEM-aero core variables

        deficits: dict(nWT) holding a list of deficits contributions from upstream wakes
        turb:  dict(nWT) holding a list of mean turbulence intensities contributions from upstream wakes
        inlets_ffor: dict(nWT) holding a list of array containing the flow field in the fixed frame of reference from
        upstream wakes contributions
        inlets_ffor_deficits: dict(nWT) holding a list of array containing the flow field in the fixed frame of
        reference from upstream wakes contributions at the rotor position
        out: dict(nWT),holding the main sDWM model outputs i.e mean power from BEM, mean power estimated from powercurve,
        mean rotor averaged wind speed, mean rotor average turbulence intensity, mean thrust coefficient from BEM and
        from power curve

        Outputs
        ----------
        aero (instance of class): updated Instance of class Aero holding BEM-aero core variables
        mfor (instance of class): updated Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
        out dict(nWT): dict including mean power from PC and BEM, mean thrust coefficient from PC and BEM
        BEM (instance of class): holds the key results from the BEM calculation

    """
    mfor   =  MFoR(meta.WTG)
    ## Compute the average wake deficit accumulation
    if meta.wake_ind_setting==1:
        rWS = deficits.get(str(meta.wtg_ind[0])) # we get the deficit of the generating plan for the current WT

        #print 'deficit of the generating plan for the current WT: ', rWS
        #raw_input('Press Enter to continue')

        if not rWS: # means that the current turbine is in free stream (first Turbine in the row)
            rWS=np.array([1.0])
        if meta.accu == 'linear':
            meta.mean_WS_DWM=meta.WS*(1.-(np.sum([1. - xx for xx in rWS])))
        elif meta.accu == 'quadratic':
            meta.mean_WS_DWM=meta.WS*(1.-np.sqrt(np.sum([(1.-xx)**2 for xx in rWS])))
        elif meta.accu == 'dominant':
            meta.mean_WS_DWM= meta.WS*(1.-(np.max([1. - xx for xx in rWS])))
            # meta.mean_WS_DWM= meta.WS*(1.-(np.max([1. - xx for xx in rWS])))
        # elif meta.accu == 'bypassed':
        #     meta.mean_WS_DWM= meta.WS
        elif meta.accu == 'ewma':
            print 'This model is currently being developped and implemented soon'
        else:
            print 'You have not specified any wake accumulation procedure... or mispelled it. Options are ''linear'', ''quadratic'' , ''dominant'' or ''ewma'''
    else:
        meta.mean_WS_DWM   = meta.WS # free stream everywhere similar to current implementation of HAWC2

    # Set buildup of turbulence
    if meta.Tbuildup_setting==1:
        ti = turb.get(str(meta.wtg_ind[0]))
        if not ti:
            meta.mean_TI_DWM  = meta.TI
        else:
            meta.mean_TI_DWM  = np.max(ti)
    else:
        meta.mean_TI_DWM  = meta.TI

    # Run BEM at accumulated deficit
    aero,BEM   =  DWM_rotor_aero(meta,aero,ID_waked)
    # domain induction
    a_domain     = np.interp(meta.vr_m,np.hstack(([aero.r_w, aero.r_w[-1]+0.01, aero.r_w[-1]+0.02])), np.hstack((( 1.0 - aero.U_w), [0., 0.])))

    if meta.BEM_AINSLIE_plot:
        print 'Aero CT: ', aero.CT
        plt.figure('Comparison between a_domain (Ainslie) and a (BEM)')
        plt.title('Comparison between a_domain (Ainslie) and a (BEM)')
        plt.plot(a_domain, meta.vr_mixl, label='a_domain')
        plt.plot(BEM.a, aero.r_t, label='a (BEM)')
        plt.xlabel('Induction Factor'), plt.ylabel('radius [R]'), plt.legend()
        plt.show()
        print ' '

    ## Compute the accumulated flow field for accurate inlet definition of the MFoR wake calculation
    # center all disks before wake accumulation
    if not inlets_ffor.get(str(meta.wtg_ind[0])):# if turbine in FREESTREAM, then mfor.Uinit is initialized in the calc_mixL module
        # mfor.U_init =   None # set to none for further initialization
        radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        mfor.U_init=f3(meta.vr_mixl)
        mfor.U_init=smooth( mfor.U_init,window_len=5)
    elif meta.accu_inlet is False:
        # mfor.U_init =   None
        radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        mfor.U_init=f3(meta.vr_mixl)
        mfor.U_init=smooth( mfor.U_init,window_len=5)
    else:   # if turbine not in the freestream, we need to compute the proper accumulated inlet to the turbine
        ranger=np.linspace(-1.,1.,meta.dR*2.)  # np.linspace(-2.,2.,meta.dR*4.)
        inlets_ffor_deficits_np_3D=np.ones((len(ranger) ,len(ranger)   , len(inlets_ffor[str(meta.wtg_ind[0])])))
        # grid_x, grid_y = np.mgrid[-1.:1.:meta.dR*2j, -1.:1.:meta.dR*2j]
        grid_x, grid_y = np.mgrid[-1.:1.:meta.dR*2j, -1.:1.:meta.dR*2j]
        for ii in range(len(inlets_ffor[str(meta.wtg_ind[0])])):
            offsets=2.*(min(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0])+abs(max(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0])-min(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0]))/2.)
            # need to interp on a new array of equal size
            values=inlets_ffor_deficits[str(meta.wtg_ind[0])][ii]
            X, Y = np.meshgrid(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0]-offsets,inlets_ffor[str(meta.wtg_ind[0])][ii][0][1])
            points=np.vstack((np.ravel(X),np.ravel(Y)))
            # wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='linear')
            #print np.shape(points), np.shape(values)
            #raw_input('in Aero, press any key to continue')
            wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='linear')
            inlets_ffor_deficits_np_3D[:,:,ii]=wake_i
            # print wake_i
        if meta.accu == 'linear':
            U_init=1.-(np.sum(1.-inlets_ffor_deficits_np_3D,axis=2))
        elif meta.accu == 'quadratic':
            U_init=1.-np.sqrt(np.sum((1.-inlets_ffor_deficits_np_3D)**2,axis=2))
        elif meta.accu == 'dominant':
            U_init=np.amin(inlets_ffor_deficits_np_3D, axis=2)
        elif meta.accu == 'ewma':
            print 'This model is currently being developped and implemented soon'
        else:
            print 'You have not specified any wake accumulation procedure... or mispelled it. Options are ''linear'', ''quadratic'' , ''dominant'' or ''ewma'''
        # Transform to axisymmetric profile inlet
        r_dist_2= np.sqrt(grid_x**2 + grid_y**2 )  #Find distance to centre of wake plane
        ffor.WS_axial_sym      = np.ones((len(np.arange(0,meta.dR+1.,1))))
        ffor.WS_axial_sym[0]=np.nanmean(U_init[r_dist_2 < (1.05*np.amin(r_dist_2))])
        #print 'meta.dR: ', meta.dR
        for i_r_pos in np.arange(1,meta.dR+1,1):
            #print 'i_r_pos: ', i_r_pos
            a=r_dist_2 > ((i_r_pos+1-1.5)*(1.0/meta.dR))# rotor neg boundaries
            bb=r_dist_2 < ((i_r_pos+1-0.5)*(1.0/meta.dR)) #rotor pos boundaries
            c=np.logical_and(a,bb)
            bin_filter = c
            tmp_ffor_flow_field_ws_mean          = U_init[bin_filter]
            ffor.WS_axial_sym[i_r_pos]           = np.nanmean(tmp_ffor_flow_field_ws_mean)
        ffor.r_sym=np.arange(0,meta.dR+1.,1)/meta.dR
        # Update the DWM inlet
        if ffor.r_sym[-1] >= meta.vr_mixl[-1]:
            # mfor.U_init = (1.0-a_domain) * np.interp(meta.vr_mixl,ffor.r_sym,ffor.WS_axial_sym)
            print 'spline interpolation'
            mfor.U_init = (1.0-a_domain) * interpolate.InterpolatedUnivariateSpline(meta.vr_mixl,ffor.r_sym,ffor.WS_axial_sym)
        else:
            print 'meta.lr_mixl: ', meta.lr_mixl
            mfor.U_init = (1.0-a_domain) * np.hstack((ffor.WS_axial_sym.ravel(), np.ones((((meta.dR * int(meta.lr_mixl))-ffor.WS_axial_sym.size),1)).ravel()))
        # Finishing
        mfor.U_init_raw=mfor.U_init
        mfor.U_init[mfor.U_init < 0.0]=0.0 # prevent from negative velocities on linear summation
        mfor.U_init=smooth( mfor.U_init,window_len=5)

    if meta.BEM_AINSLIE_plot:
        print 'U_init shape: ', np.shape(mfor.U_init)
        plt.figure('U_init, the Input for Ainslie Coming from BEM')
        plt.title('U_init, the Input for Ainslie Coming from BEM')
        plt.plot(mfor.U_init, meta.vr_mixl, label='U_init')
        plt.xlabel('r [R]'), plt.ylabel('U[U0]'), plt.legend()
        plt.show()
        print ' '

    # Power curve based
    try:
        if BEM.derated is False:
            # print 'use standard ws for curve'
            aero.pow_cur=meta.WTG_spec.get_P(meta.mean_WS_DWM)
            aero.ct_cur=meta.WTG_spec.get_CT(meta.mean_WS_DWM)
        else:
            # print 'use demanded ws for curve'
            aero.pow_cur=meta.WTG_spec.get_P(BEM.Ud)
            aero.ct_cur=meta.WTG_spec.get_CT(BEM.Ud)
    except:
        aero.pow_cur=0.
        aero.ct_cur=0.
    # write outlets
    out[str(meta.wtg_ind[0])]=[]
    out[str(meta.wtg_ind[0])].append(float(format(aero.Power/1000., '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(meta.mean_WS_DWM, '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(meta.mean_TI_DWM, '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(aero.CT/1., '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(aero.pow_cur, '.2f'))) # based on power curve
    out[str(meta.wtg_ind[0])].append(float(format(aero.ct_cur, '.2f'))) # based on power curve
    return aero, mfor, out, BEM


def DWM_rotor_aero(meta,aero,ID_waked, *args):
    """
    Function that calls the BEM calculation and calculate the corrected initial wake radius
    (boundary condition to Ainslie model ref [2])

    Parameters
    ----------
    aero    class holding the aero parameters to the Ainslie model
    meta    class holding grid and ambient parameters

    Returns
    -------
    aero (updated)
             CP [1],        : turbine power coefficient
             CPloc[r,1],    : local CP for power estimation
             Power [1],     : power from BEM [W]
             CT [1],        : thrust coefficient [-]
             RPM [1],       : rotational speed [RPM]
             PITCH [1]      : pitch angle [deg]
             U_w [1,r],     : local axial velocity in the wake from momemtum theory
             a [1,r],       : local induction factor
             dA [1],        : local annular area
             f_w [1],       : calibration factor on wake deficit expansion
             mean_a [1],    : rotor averaged axial induction
             r_t [1,r],     : non dimensional radial position where a is evaluated
             r_w [1,r]      : wake radius in near wake regime
    """

    if (meta.mean_WS_DWM >= meta.WTG_spec.u_cutin) or (meta.mean_WS_DWM <= meta.WTG_spec.u_cutout) is True:
        BEM =getInduction(30, meta.WTG, 'hawc', meta.mean_WS_DWM, meta,ID_waked, derating=meta.derating)
        aero.a = np.array(BEM.a)
        aero.r_t = np.array(BEM.r)/BEM.R
        aero.CP = np.array(BEM.CP)
        # aero.CPloc = np.array(BEM.CPloc)
        aero.Power = np.array(BEM.Power)
        aero.CT    = np.array(BEM.CT)

        aero.RPM    = BEM.RPM
        aero.RPM_opt   = BEM.RPM_opt
        aero.PITCH    = BEM.PITCH
        # print 'pitch is',aero.PITCH
        aero.PITCH_opt    = BEM.PITCH_opt
        # aero.a = np.hstack((0.,aero.a,0.))
        # aero.r_t = np.hstack((0.,aero.r_t,1.))
    else:
        BEM = InitBEM(30)
        aero.r_t = np.arange(0.,1.+1./(meta.dR),1./(meta.dR))
        aero.Power  = 0.0 #kW
        aero.CT     = meta.WTG_spec.CT_idle
        aero.ia = meta.WTG_spec.get_a(aero.CT)
        aero.a=  aero.ia*np.ones(len(aero.r_t))
    # print 'BEM predicts %4.2f kW at %4.2f m/s' %(BEM.Power/1000.,meta.mean_WS_DWM)

    # Boundary conditions for the r_w
    aero.dA = np.concatenate(([0], (pi*aero.r_t [1:]**2 - pi*aero.r_t [0:-1]**2)), axis=0)
    aero.mean_a= np.sum(aero.a*aero.dA)/pi
    # Uniform expansion
    aero.f_w        = sqrt( (1.0-aero.mean_a) / (1.0- ((1.0+meta.fR) * aero.mean_a))  )
    aero.r_w = np.dot(aero.r_t, aero.f_w)
    # Boundary conditions for U_w
    aero.U_w     = 1.0-(aero.a * (1.0 + meta.fU))

    return aero, BEM


def DWM_make_inflow_to_mixl(meta,ffor,inlets):
    """
    Function that performs the velocity integration for each downstream rotor of a given upstream turbine.

    Parameters
    ----------
    ffor:(instance of class): Instance of class Ffor holding the fixed frame of reference velocity field in global
    WF coordinates
    meta (instance of class): Instance of class Meta holding DWM core variables
    deficits: dict(nWT) holding a list of deficits contributions from upstream wakes
    ID_waked: dict(nWT) holding list of upstream turbine index for each turbine in the wind farm
    inlets_ffor: dict(nWT) holding a list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_deficits: dict(nWT) holding a list of array containing the flow field in the fixed frame of reference from upstream wakes contributions at the rotor position

    Returns
    -------
    deficits: dict(nWT) updated list of deficits contributions from upstream wakes
    ID_waked: dict(nWT) updated list of upstream turbine index for each turbine in the wind farm
    inlets_ffor: dict(nWT) updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_deficits: dict(nWT) updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions at the rotor position
    """

    # meta.nr                = round(np.sqrt((max(meta.x_vec - meta.hub_x[0]))**2 + (max(meta.y_vec - meta.hub_y))**2 ) * meta.dR)
    meta.nr   =  max([round(np.sqrt((max(meta.x_vec - meta.hub_x[i_z]))**2 + (max(meta.y_vec - meta.hub_y))**2 ) * meta.dR) for i_z in np.arange(0,meta.nz,1)])
    max_index = [round(np.sqrt((max(meta.x_vec - meta.hub_x[i_z]))**2 + (max(meta.y_vec - meta.hub_y))**2 ) * meta.dR) for i_z in np.arange(0,meta.nz,1)].index(meta.nr)

    ffor.WS_axial_sym      = np.ones((meta.nr,meta.nz))
    # ffor.WS_axial_sym_pow3 = np.ones((meta.nr,meta.nz))
    ffor.TI_meand_axial_sym = np.ones((meta.nr,meta.nz))
    ffor.TI_axial_sym      = np.ones((meta.nr,meta.nz))
    for i_z in np.arange(0,meta.nz,1):
        meta.r_dist_2= np.sqrt((meta.x_mat - meta.hub_x[i_z])**2 + (meta.y_mat - meta.hub_y)**2 )  #Find distance to centre of wake plane
        tmp_ws       = np.squeeze(ffor.WS_axial_ffor[:,:,i_z])
        # tmp_ws_pow3  = np.squeeze(ffor.WS_axial_ffor_pow3[:,:,i_z])
        tmp_meand_TI = np.squeeze(ffor.TI_meand_axial_ffor[:,:,i_z])
        tmp_TI       = np.squeeze(ffor.TI_axial_ffor[:,:,i_z])
        ffor.WS_axial_sym[0,i_z]       = np.mean(tmp_ws[meta.r_dist_2 < (0.5)*(1.0/meta.dR)])
        # ffor.WS_axial_sym_pow3[0,i_z]  = np.mean(tmp_ws_pow3[meta.r_dist_2 < (0.5)*(1.0/meta.dR)])
        ffor.TI_axial_sym[0,i_z]       = np.mean(tmp_TI[meta.r_dist_2 < (0.5)*(1.0/meta.dR)])
        ffor.TI_meand_axial_sym[0,i_z]  = np.mean(tmp_meand_TI[meta.r_dist_2 < (0.5)*(1.0/meta.dR)])

        for i_r_pos in np.arange(1,round(np.sqrt((max(meta.x_vec - meta.hub_x[i_z]))**2 + (max(meta.y_vec - meta.hub_y))**2 ) * meta.dR),1):
            a=meta.r_dist_2 > ((i_r_pos+1-1.5)*(1.0/meta.dR))# rotor neg boundaries
            bb=meta.r_dist_2 < ((i_r_pos+1-0.5)*(1.0/meta.dR)) #rotor pos boundaries
            c=np.logical_and(a,bb)
            bin_filter = c
            tmp_ffor_flow_field_ws_mean          = tmp_ws[bin_filter]
            # tmp_ffor_flow_field_ws_mean_pow3     = tmp_ws_pow3[bin_filter]
            tmp_ffor_flow_field_TI_mean          = tmp_TI[bin_filter]
            tmp_ffor_flow_field_meand_TI_mean    = tmp_meand_TI[bin_filter]
            ffor.WS_axial_sym[i_r_pos,i_z]     = np.mean(tmp_ffor_flow_field_ws_mean)
            # ffor.WS_axial_sym_pow3[i_r_pos,i_z] = np.mean(tmp_ffor_flow_field_ws_mean_pow3)
            ffor.TI_axial_sym[i_r_pos,i_z]      = np.mean(tmp_ffor_flow_field_TI_mean**2)**0.5
            ffor.TI_meand_axial_sym[i_r_pos,i_z] = np.mean(tmp_ffor_flow_field_meand_TI_mean**2)**0.5

        ffor.r_sym=np.arange(0,round(np.sqrt((max(meta.x_vec - meta.hub_x[max_index]))**2 + (max(meta.y_vec - meta.hub_y))**2 ) * meta.dR),1)/meta.dR
        inlets[str(meta.wtg_ind[i_z])].append([np.vstack((ffor.r_sym,ffor.WS_axial_sym[:,i_z]))])

        # inlets_ffor[str(meta.wtg_ind[i_z])].append([np.vstack(((meta.x_vec-meta.hub_x[i_z])/2.,ffor.WS_axial_ffor[:,:,i_z]))])
        # if str(meta.wtg_ind[i_z])=='18' and str(meta.wtg_ind[0])=='18':
        #     print 'inlets', [np.vstack((ffor.r_sym,ffor.WS_axial_sym[:,i_z]))]

    return ffor,inlets


def DWM_get_deficit(ffor,meta,deficits,ID_waked,inlets_ffor,inlets_ffor_deficits):

    #raw_input('Begin to get deficit')
    for i_z in np.arange(0,meta.nz,1):
        # on global frame mesh
        X,Y = np.meshgrid(ffor.x_vec, ffor.y_vec)
        index_trapz=np.sqrt((X + meta.C2C[i_z]/(2.*meta.WTG_spec.R))**2 + (Y)**2 )>=0.5
        wakedefmask = np.ma.array(np.squeeze(ffor.WS_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=0.0).filled()
        wakedefmasknancoarse = np.ma.array(np.squeeze(ffor.WS_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=np.nan).filled()
        disk = np.ma.array(np.zeros(wakedefmask.shape), mask=~index_trapz, fill_value=1.0).filled()
        disk_area=np.trapz(np.trapz(disk,dx=1./meta.dy),dx=1./meta.dx)
        trapz2=np.trapz(np.trapz(wakedefmask,dx=1./meta.dy),dx=1./meta.dx)


        # on finer mesh
        # values=np.squeeze(ffor.WS_axial_ffor[:,:,i_z])
        # X, Y = np.meshgrid(ffor.x_vec,ffor.y_vec)
        # points=np.vstack((np.ravel(X),np.ravel(Y)))
        # grid_x, grid_y = np.mgrid[-0.5:0.5:meta.dR*36j, -0.5:0.5:meta.dR*36j]
        # wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='cubic')
        #
        # X2,Y2=np.meshgrid(np.linspace(-0.5,0.5,36*meta.dR),np.linspace(-0.5,0.5,36*meta.dR))
        # index_trapz=np.sqrt((X2 + meta.C2C[i_z]/(2.*meta.WTG_spec.R))**2 + (Y2)**2 )>=0.5
        #
        # wakedefmask = np.ma.array(wake_i, mask=index_trapz, fill_value=0.0).filled()
        # # wakedefmasknan = np.ma.array(wake_i, mask=index_trapz, fill_value=np.nan).filled()
        # disk = np.ma.array(np.zeros(wakedefmask.shape), mask=~index_trapz, fill_value=1.0).filled()
        # disk_area=simps(simps(disk,np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dy)),np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dx))
        # trapz2=simps(simps(wakedefmask,np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dy)),np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dx))
        #
        #print wakedefmasknancoarse

        #Updating
        deficits[str(meta.wtg_ind[i_z])].append(trapz2/disk_area)
        ID_waked[str(meta.wtg_ind[i_z])].append(meta.wtg_ind[0])
        inlets_ffor_deficits[str(meta.wtg_ind[i_z])].append(wakedefmasknancoarse)
        inlets_ffor[str(meta.wtg_ind[i_z])].append([np.vstack(((meta.x_vec-meta.hub_x[i_z]),(meta.y_vec-meta.hub_y),ffor.WS_axial_ffor[:,:,i_z]))])
        #print 'inlets_ffor shape (after mean): ', np.shape([np.vstack(((meta.x_vec-meta.hub_x[i_z]),(meta.y_vec-meta.hub_y),ffor.WS_axial_ffor[:,:,i_z]))])
        #print 'inlets_ffor_deficits shape (after mean): ', np.shape(wakedefmasknancoarse)

    #for i_z in range(0, meta.nz):
        #print 'Deficit for WT i_z = ', str(i_z), ': ', deficits[str(i_z)]
        #print 'inlets ffor Deficit for WT i_z = ', str(i_z), ': ', inlets_ffor_deficits[str(meta.wtg_ind[i_z])]
        #print 'inlets for WT i_z = ', str(i_z), ': ', inlets_ffor[str(meta.wtg_ind[i_z])]


    # Plotting
    # Useless to plot Inlet, it comes directly from WS_axial plotted in meandering part plus some coordinate vector
    """
    plt.figure('Inlets FFoR')
    length_ref = (len(deficits[str(0)])-1)
    plt.title('Inlet generated by the WT'+str(length_ref))
    #plt.contourf(inlets_ffor[str(meta.wtg_ind[0])][length_ref][0])
    plt.pcolor(inlets_ffor[str(meta.wtg_ind[0])][length_ref][0])
    plt.colorbar()
    plt.show()
    #"""
    if meta.DEFICIT_plot:
        plt.figure('Deficit')
        plt.title('Deficits generating by each Turbine on other Turbines'), plt.xlabel('WT'), plt.ylabel('Deficit')
        for i in range(len(deficits[str(0)])):
            length_ref = (len(deficits[str(0)])-1) + meta.nz - i
            #print 'i=', i
            #print 'length ref ', length_ref
            Deficit_to_plot = [deficits[str(i_z)][i] for i_z in np.arange(0, length_ref, 1)]
            #print Deficit_to_plot
            plt.plot(np.arange(0, length_ref, 1), Deficit_to_plot, label='WT'+str(length_ref-1))
            if i==0:
                plt.xlim(length_ref-1, 0)
        plt.legend()
    #raw_input('End of Get_deficit')
    return deficits, ID_waked,inlets_ffor,inlets_ffor_deficits


def DWM_get_turb(ffor,meta,turb,inlets_ffor_turb,):
    """
    Function that calculate the rotor averaged turbulence intensity
    Parameters
    ----------
    ffor:(instance of class): Instance of class Ffor holding the fixed frame of reference velocity field in global
    WF coordinates
    meta (instance of class): Instance of class Meta holding DWM core variables
    turb:  dict(nWT) holding a list of turbulence intensities contributions from upstream wakes
    inlets_ffor_turb: dict(nWT) holding a list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position

    Returns
    -------
    turb:  dict(nWT) updated list of turbulence intensities contributions from upstream wakes
    inlets_ffor: dict(nWT) updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_turb: dict(nWT) updated list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position
    """
    for i_z in np.arange(0,meta.nz,1):
        X, Y = np.meshgrid(ffor.x_vec, ffor.y_vec)
        index_trapz=np.sqrt((X + meta.C2C[i_z]/(2.*meta.WTG_spec.R))**2 + (Y)**2 )>=0.5
        turbmask = np.ma.array(np.squeeze(ffor.TI_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=0.0).filled()
        turbmasknan = np.ma.array(np.squeeze(ffor.TI_axial_ffor[:,:,i_z]), mask=index_trapz, fill_value=np.nan).filled()
        disk = np.ma.array(np.zeros(turbmask.shape), mask=~index_trapz, fill_value=1.0).filled()
        disk_area=np.trapz(np.trapz(disk,dx=1./meta.dy),dx=1./meta.dx)
        trapz2=np.trapz(np.trapz(turbmask,dx=1./meta.dy),dx=1./meta.dx)
        turb[str(meta.wtg_ind[i_z])].append(trapz2/disk_area)
        inlets_ffor_turb[str(meta.wtg_ind[i_z])].append(turbmasknan)
        #print 'inlets_ffor_deficits shape (after mean): ', np.shape(turbmasknan)

    if meta.DEFICIT_plot:
        plt.figure('Rotor averaged turbulence intensity')
        plt.title('Rotor averaged turbulence intensity generating by each Turbine on other Turbines'), plt.xlabel('WT'), plt.ylabel('Deficit')
        for i in range(len(turb[str(0)])):
            length_ref = (len(turb[str(0)])-1) + meta.nz - i
            #print 'i=', i
            #print 'length ref ', length_ref
            Turb_to_plot = [turb[str(i_z)][i] for i_z in np.arange(0, length_ref, 1)]
            #print Deficit_to_plot
            plt.plot(np.arange(0, length_ref, 1), Turb_to_plot, label='WT'+str(length_ref-1))
            if i==0:
                plt.xlim(length_ref-1, 0)
        plt.legend()
        plt.show()
    return turb,inlets_ffor_turb


def DWM_MFOR_to_FFOR(mfor,meta,meand,ffor):
    """
    Function that calculate the velocity in the fixed (global) frame of reference from the Mfor

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
    ffor.ffor_flow_field_TI_tmp_tmp =  meta.TI * np.ones((meta.nx, meta.ny))  #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.TI_axial_ffor_tmp     =  np.zeros((meta.nx, meta.ny, meta.nz))        #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.WS_axial_ffor_tmp     =  np.zeros((meta.nx, meta.ny, meta.nz))        #X = lateral ,Y = vertical, time ,Z = streamwise
    ffor.ffor_flow_field_ws_tmp2    =  np.zeros((meta.nx, meta.ny, meta.nz))   #X = lateral ,Y = vertical, time ,Z = streamwise
    #print 'shape(ffor.WS_axial_ffor_tmp): ', np.shape(ffor.WS_axial_ffor_tmp)
    ffor.TI_meand_axial_ffor    =  np.zeros((meta.nx, meta.ny, meta.nz))
    ffor.WS_axial_ffor =  np.zeros((meta.nx, meta.ny, meta.nz))
    ffor.TI_axial_ffor =  np.zeros((meta.nx, meta.ny, meta.nz))
    #print 'shape(ffor.WS_axial_ffor): ', np.shape(ffor.WS_axial_ffor)
    ffor.x_vec_t            = np.zeros((meta.nx,meta.nz))
    ffor.x_mat_t            = np.zeros((meta.nx,meta.ny,meta.nz))
    # Define radial distance vectors
    # r_dist_2              = 999*np.ones((meta.nx,meta.ny))
    r_dist                = 999*np.ones((meta.nx,meta.ny))

    # CREATES THE GLOBAL FLOW FIELDS IN CARTESIAN GRID
    # meta.x_mat = np.tile(meta.x_vec.reshape(len(meta.x_vec),1),meta.ny)
    # meta.y_mat = np.tile(meta.y_vec,(meta.nx,1))
    meta.z_mat = np.tile(meta.vz,(meta.nx,1))
    #print 'meta vz: ', meta.vz

    # Store the ffor flow field
    ffor.x_vec = (meta.x_vec - meta.hub_x[0]) / 2.
    ffor.y_vec = (meta.y_vec - meta.hub_y) / 2.
    ffor.z_vec = meta.z_vec + np.hstack((0., np.cumsum(meta.hub_z[0:])))[0]
    ffor.x_mat = meta.x_mat / 2.        # [R] -> [D]
    ffor.y_mat = meta.y_mat / 2.        # [R] -> [D]
    ffor.z_mat = meta.z_mat / meta.dz   # [D]

    if meta.MEANDERING_plot:
        plt.ion()
        plt.figure()

    for i_z in np.arange(0, meta.nz, 1):

        #print 'meta.vz[iz]: ', meta.vz[i_z]
        # EXTRACT TI_DWM AND WS_DWM IN MFoR
        try:
            DWM_WS_DATA = mfor.U[meta.vz[i_z], :]

            #print 'DWM_WS_DATA: ', DWM_WS_DATA


        except:
            print 'Fatal error, possibly due to a too low turbulence intensity with respect to the demanded mean wind speed, try increase the input TI'
        DWM_TI_DATA = mfor.TI_DWM[meta.vz[i_z],:]
        #print 'DWM_TI_DATA: ', DWM_TI_DATA

        ### Correct DWM_TI_DATA so that no point to have lower TI than "TIamb"
        DWM_TI_DATA[DWM_TI_DATA < np.nanmean(meta.mean_TI_DWM)] = np.nanmean(meta.mean_TI_DWM)





        #for i_t in np.arange(0, 2, 1):
        for i_t in np.arange(0,len(meand.time),1):
            Ro_x                    = meand.meand_pos_x[i_z,i_t]
            Ro_y                    = meand.meand_pos_y[i_z,i_t]

            #print '(Ro_x,Ro_y): ', [Ro_x, Ro_y]
            #print 'meta.x_mat: ', meta.x_mat
            r_dist                  = np.sqrt((meta.x_mat - Ro_x)**2 + (meta.y_mat - Ro_y)**2 )
            #print 'r_dist: ', r_dist
            #print 'r_dist: ', r_dist

            ###############################################################
            #do we have to keep this for dynamic, or use GLarsen function?
            #print 'mfor.wakeW[meta.vz[i_z]]: ', mfor.WakeW[meta.vz[i_z]]
            tmp_index               = r_dist < mfor.WakeW[meta.vz[i_z]]*1.5
            ################################################################
            #print 'tmp_index: ', tmp_index
            #raw_input('entry')
            tmp_field_WS            = np.ones((meta.nx,meta.ny))
            #print 'tmp_index: ', tmp_index

            # it's here that we change the velocity to be in FFOR
            tmp_field_WS[tmp_index] = np.interp(r_dist[tmp_index], meta.vr_m, DWM_WS_DATA)


            ffor.WS_axial_ffor_tmp[:, :, i_z]  = ffor.WS_axial_ffor_tmp[:, :, i_z] + (tmp_field_WS)
            ffor.ffor_flow_field_ws_tmp2[:, :, i_z] = ffor.ffor_flow_field_ws_tmp2[:, :,i_z] + (tmp_field_WS**2)

            tmp_field_TI            = meta.TI * np.ones((meta.nx,meta.ny))
            tmp_field_TI[tmp_index] = np.interp( r_dist[tmp_index],meta.vr_m,DWM_TI_DATA)

            ffor.ffor_flow_field_TI_tmp_tmp[:, :]      = tmp_field_TI
            ffor.TI_axial_ffor_tmp[:, :, i_z]     = ffor.TI_axial_ffor_tmp[:, :, i_z] + ffor.ffor_flow_field_TI_tmp_tmp**2

            if meta.MEANDERING_plot:
                plt.subplot(131)
                plt.title('Axial Velocity at WT ' + str(7 - i_z)+' in MFoR')
                plt.plot(DWM_WS_DATA, meta.vr_mixl, label='DWM_WS_DATA')
                plt.xlabel('U (MFoR)'), plt.ylabel('r [R]')

                plt.subplot(132)
                plt.title('Meandering WS in Time, at WT ' + str(7 - i_z) + ' in FFoR')
                CF = plt.contourf(ffor.x_mat, ffor.y_mat, tmp_field_WS, np.arange(0.2, 1, .05), extend='both')
                plt.xlabel('Lateral direction, x [D]'), plt.ylabel('Longitudinal direction, y [D]')
                plt.colorbar(CF)

                plt.subplot(133)
                plt.title('Meandering TI in Time, at WT ' + str(7 - i_z) + ' in FFoR')
                CF = plt.contourf(ffor.x_mat, ffor.y_mat, tmp_field_TI, np.arange(0.08, 0.4, .02), extend='both')
                plt.xlabel('Lateral direction, x [D]'), plt.ylabel('Longitudinal direction, y [D]')
                plt.colorbar(CF)

                plt.draw()
                plt.pause(0.1)
                plt.cla()
                plt.clf()
        #plt.plot([0,100],[0,1],ffor.WS_axial_ffor_tmp[:, :, i_z])
        # """

    if meta.MEANDERING_plot:
        plt.ioff()

        #"""
        #"""
    # Stores the mean field
    for i_z in np.arange(0,meta.nz,1):
        #### Here we average all time related data along the z-axis
        # don't keep this part for dynamic
        ffor.TI_meand_axial_ffor[:, :, i_z]=np.sqrt(abs(ffor.ffor_flow_field_ws_tmp2[:, :, i_z] - ((ffor.WS_axial_ffor_tmp[:, :, i_z]**2)/len(meand.time)) )/ (len(meand.time)-1.0))
        ffor.WS_axial_ffor[:, :, i_z]      = (ffor.WS_axial_ffor_tmp[:, :, i_z]  / len(meand.time))
        ffor.TI_axial_ffor[:, :, i_z]      = (ffor.TI_axial_ffor_tmp[:, :, i_z]  / len(meand.time))**(1.0/2.0)

        ffor.x_vec_t[:, i_z]               = (meta.x_vec-meta.hub_x[i_z])/2.
        ffor.x_mat_t[:, :, i_z]              = np.tile(ffor.x_vec_t[:, i_z] .reshape(len(ffor.x_vec_t[:, i_z]),1),meta.ny)/2.

    if meta.MEANDERING_plot:
        for i_z in np.arange(0, meta.nz, 1):

            plt.figure('Averaged Meandering WS for statistical approach (FFoR) at WT'+str(7-i_z))
            plt.title('Averaged axial WS Field at WT' + str(7 - i_z))
            CF1 = plt.contourf(ffor.x_mat, ffor.y_mat, ffor.WS_axial_ffor[:, :, i_z])
            plt.xlabel('Lateral direction, x [D]'), plt.ylabel('Longitudinal direction, y [D]')
            plt.colorbar(CF1)


            plt.figure('Averaged Meandering TI for statistical approach (FFoR) at WT' + str(7 - i_z))
            plt.subplot(121)
            plt.title('Averaged axial TI at WT' + str(7 - i_z))
            CF2 = plt.contourf(ffor.x_mat, ffor.y_mat, ffor.TI_axial_ffor[:, :, i_z])
            plt.xlabel('Lateral direction, x [D]'), plt.ylabel('Longitudinal direction, y [D]')
            plt.colorbar(CF2)

            plt.subplot(122)
            plt.title('Averaged axial meandering TI at WT' + str(7 - i_z))
            CF3 = plt.contourf(ffor.x_mat, ffor.y_mat, ffor.TI_meand_axial_ffor[:, :, i_z])
            plt.xlabel('Lateral direction, x [D]'), plt.ylabel('Longitudinal direction, y [D]')
            plt.colorbar(CF3)
        plt.show()

    return mfor,ffor,meta,meand


def DWM_make_grid(meta):
    """
    Function to adapt the grid (polar MFoR and cartesian FFoR) length
    in the streamwise direction based on the distance to the next turbine
    It also create extraction plane in the cartesian ffor domain

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
    #print 'meta.vz in makegrid(float): ', meta.vz
    meta.vz = np.rint(meta.vz).astype(dtype=int)
    #print 'meta.vz rounded: ', meta.vz
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


def DWM_meta_meand(meand,meta):
    """
    Function to calculate the standard deviation of the wake center in-plane.

    Parameters
    ----------
    meand   class holding the meandering parameters
    meta    class holding grid and ambient parameters

    Returns
    -------
    meand (updated)
            meand_pos_x [z,t]  : vector holding wake position in lateral
            meand_pos_y [z,t]  : vector holding wake position in longitunal
            std_meand_x [z,t]  : standard deviation of wake center position in lateral
            std_meand_y [z,t]  : standard deviation of wake center position in long
            time   [t]   : time vector
    """
    #meand.time             = np.arange(1,200.+1,1) # number of wake meandering samples 10 min at 1Hz
    meand.time             = np.arange(1,100+1,1) # number of wake meandering samples 10 min at 1Hz
    # meand.time         = np.arange(1,500+1,1) # number of wake meandering samples 10 min at 1Hz for DEBUG ONLY
    # meand.x_offset     = 2.0*np.tan(meta.dir*(pi/180.)) * meta.z_vec #offset to wind dir relative row

    # Wake meandering based on meandering method with 0.8 as wake transport speed, see Keck et al. [4] & [5]
    meand.std_meand_x, meand.std_meand_y = meand_table_DWM_method(meta)

    # build meandering vectors at the planes specified by meta.vz
    meand.meand_pos_x=np.zeros((meta.nz,len(meand.time)))
    meand.meand_pos_y=np.zeros((meta.nz,len(meand.time)))

    # FOR DEBUG ONLY
    # meand.rea = io.loadmat('realizations_debug/realization.mat')
    # here the meandering patterns are assumed independent from each other between downstream planes
    # this is not physical as meandering paths are well correlated between consecutives downstream distances
    # with a scaling = variation of std_meand_x
    # time shift due to advection time
    # seed_x=np.random.randn(len(meand.time),1)
    # seed_y=np.random.randn(len(meand.time),1)
    if meta.optim is True:
        # meand.rea = io.loadmat('realizations_debug/realization.mat')
        # seed_x=meand.rea.get('realization')
        # seed_y=meand.rea.get('realization')
        seed_x=np.load('../data/meand_x.npy')
        seed_y=np.load('../data/meand_y.npy')


        #print 'seed_x: ', seed_x
        #print 'len(seed_x): ', len(seed_x)
        #print 'seed_y: ', seed_y
        #print 'len(seed_y): ', len(seed_y)
    elif meta.optim is False:
        seed_x=np.random.randn(len(meand.time),1)
        seed_y=np.random.randn(len(meand.time),1)
    #seed_x=meand.rea.get('realization')
    #seed_y=meand.rea.get('realization')
    #####PLOT Meandering####
    #print seed_x
    #print meta.hub_x[0]
    #print meand.std_meand_x
    #print meand.meand_pos_x[1,0]
    if meta.MEANDERING_plot:
        for i_z in np.arange(meta.nz-1, -1, -1):
            plt.figure('Wake Center position (statistical approach) for WT '+str(7-i_z))
            plt.title('Wake Center position (statistical approach) for WT '+str(7-i_z))
            plt.plot(meand.meand_pos_x[i_z,:], meand.meand_pos_y[i_z,:], 'x')
            plt.xlabel('Lateral deviation, x [R]'), plt.ylabel('Longitudinal deviation, y [R]')
            plt.xlim((-6, 9))
            plt.ylim((-2, 8))
        plt.show()
    for i_z in np.arange(0, meta.nz, 1):
        meand.meand_pos_x[i_z, :] = (meta.hub_x[0] + (meand.std_meand_x[i_z] * seed_x)).ravel()
        meand.meand_pos_y[i_z, :] = (meta.hub_y + (meand.std_meand_y[i_z] * seed_y)).ravel()
    return meand


def meand_table_DWM_method(meta):
    """
    Function to determine the meandering magnitude as function of ambient conditions at given downstream position i_z
    The parameters required for the determination of the meandering magnitude are: downstream distance, height,
    turbulence intensity and atmospheric stability
    Parameters
    ----------
    meta    class holding grid and ambient parameters

    Returns
    -------
    std_meand_x[i_z,1]   standard deviation of lateral meandering magnitude
    std_meand_y[i_z,1]   standard deviation of longitudinal meandering magnitude

    """
    tmp_TI=np.zeros((1,2)).ravel();tmp_iTI=np.zeros((1,2)).ravel()
    print 'location = meta.z_vec:' + str(meta.z_vec)+'  [D]'
    location=meta.z_vec
    lloc=len(location)
    print 'number of remaining wind Turbine (lloc): ', lloc   #number of remaining wind Turbine
    index_orig=np.argsort(location)
    location=np.sort(location) # sorting for interpolation on continuous function
    Meand = io.loadmat('../data/meand_data.mat')
    #print 'Meand: ', Meand

    # TI reference vector for Meand matrix
    TI_vector     = np.array([0.0, 60.0, 100.0, 140.0, 200.0])/1000.0
    #creates tmp_TI for interpolation
    ind=TI_vector<=meta.TI
    r = np.array(range(len(ind)))
    tmp_iTI[0] = max(r[ind])
    tmp_TI[0]  = max(TI_vector[r[ind]])
    try:
        ind=TI_vector>=meta.TI
        r = np.array(range(len(ind)))
        tmp_iTI[1]= min(r[ind])
        tmp_TI[1]  = min(TI_vector[r[ind]])
    except:
        tmp_TI[1]  = tmp_TI[0]
    # Test if TI over TI table limit
    if (tmp_TI[0] == tmp_TI[1])==1:
        tmp_TI_length = 1
    else:
        tmp_TI_length = 2

    #print 'tmp_TI: ', tmp_TI
    #print 'tmp_TI_length', tmp_TI_length

    std_hor_tmp=np.zeros((tmp_TI_length,lloc));std_vert_tmp=np.zeros((tmp_TI_length,lloc))
    std_meand_x=np.zeros((lloc));std_meand_y=np.zeros((lloc))
    #std_hor_tmp=[];std_vert_tmp=[]
    for i_iTI in np.arange(0,tmp_TI_length,1,'int16'):
        print 'Performing for Penas Mann stability'
        # Atmospheric stability effects based on Penas Mann stability & MARDM integration
        #print 'meta.atmo_stab', meta.atmo_stab
        #print 'i_iTI: ', i_iTI
        #print 'tmp_iTI[i_iTI]): ', tmp_iTI[i_iTI]
        #According to 2013-01-11 mail
        tmp_data_hor= Meand[meta.atmo_stab][0][0][0][0][0][ int(tmp_iTI[i_iTI])]
        tmp_data_vert= Meand[meta.atmo_stab][0][0][1][0][0][ int(tmp_iTI[i_iTI])]
        #print 'tmp data hor:', tmp_data_hor
        #print 'tmp data vert:', tmp_data_vert

        #  Adds value at D=0 and 108D (=18D*6)
        tmp_data_hor  = np.vstack(([0.0, 0.0, 0.0],tmp_data_hor,tmp_data_hor[-1,:]*6.))
        tmp_data_vert = np.vstack(([0.0, 0.0, 0.0],tmp_data_vert,tmp_data_vert[-1,:]*6.))
        #print 'tmp data hor after Adding value at D=0 and 108D (=18D*6):', tmp_data_hor
        #print 'tmp data vert after Adding value at D=0 and 108D (=18D*6):', tmp_data_vert

        #specific the format of the meandering data
        dist_vector   = np.array([0.0, 1.0, 2.0, 3.0, 4.5, 6.0, 7.5, 9.0, 12.0, 18.0, 108.0])
        height_vector = np.array([40.0, 100.0, 160.0])

        #finds the wake meandering
        sp = interpolate.RectBivariateSpline(dist_vector, height_vector, tmp_data_hor, kx=1, ky=1, s=0)
        std_hor_tmp[i_iTI,:]=sp(location,meta.WTG_spec.H).reshape(1,lloc)

        sp = interpolate.RectBivariateSpline(dist_vector, height_vector, tmp_data_vert, kx=1, ky=1, s=0)
        std_vert_tmp[i_iTI,:]=sp(location,meta.WTG_spec.H).reshape(1,lloc)

    if (tmp_TI_length == 1)==1 and (tmp_TI[1] == meta.TI)==1: # use tmp_TI(1) value
        std_meand_x = std_hor_tmp[0,:]
        std_meand_y = std_vert_tmp[0,:]
    elif (tmp_TI_length == 1)==1: # Scale tmp_TI(1) value to TI_AMB value
        std_meand_x  = (meta.TI/tmp_TI[0])* std_hor_tmp[0,:]
        std_meand_y = (meta.TI/tmp_TI[0])* std_vert_tmp[0,:]
    else: # interpolate between tmp_TI(1) and tmp_TI(2).
        for i_z in np.arange(0,lloc,1):
            std_meand_x[i_z]  = np.interp(meta.TI,[tmp_TI[0], tmp_TI[1]],std_hor_tmp[:,i_z])
            #np.interp(meta.WTG_spec.H,[40,  100,  160], L_ABL_vector  )
            std_meand_y[i_z] = np.interp(meta.TI,[tmp_TI[0], tmp_TI[1]],std_vert_tmp[:,i_z])


    std_meand_x=std_meand_x[index_orig] # reorder to original sorting
    std_meand_y=std_meand_y[index_orig]
    #print 'Std_meand_y: ', std_meand_y
    #print 'Std_meand_x: ', std_meand_x

    if meta.MEANDERING_plot:
        plt.figure('Standard deviation (statistical approach)')
        plt.title('Standard deviation (statistical approach) for each WindTurbines')
        plt.plot(range(meta.nz-1, -1, -1), std_meand_x, 'o', label= 'Std x')
        plt.plot(range(meta.nz-1, -1, -1), std_meand_y, 'x', label='Std y')
        plt.xlabel('Downstreamwise'), plt.ylabel('Standard Deviation')
        plt.show()

    return std_meand_x, std_meand_y


def DWM_outputs(DWM,ffor,mfor,meta, aero,par, BEM):
    """
    Function that store chosen flow field and turbine results

    Parameters
    ----------
    DWM  [nWT]      list container for all results
    meta              class instance holding grid and ambient parameters
    ffor              class instance holding the fixed frame of reference flow field data
    mfor              class instance holding the meandering frame of reference flow field data
    Returns
    -------
    DWM (updated)
       For the description of the returned variables, refer to the class definition.


             WS_axial_ffor [x,y,z]        : streamwise wake velocity deficit [-]
             WS_axial_ffor_pow3 [x,y,z]   : streamwise wake velocity deficit ^3 for power estimation Eq [14] [-]
             TI_axial_ffor    [x,y,z]     : apparent turbulence intensity in FFoR [-]
             TI_meand_axial_ffor [x,y,z]  : turbulence intensity contribution of wake meandering Eq [10] [-]
             TI_tot_axial_ffor [r,z]      : axi-symmetric total turbulence intensity [-]
             WS_axial_sym [r,z]           : axi-symmetric radial velocity deficit [-]
             WS_axial_sym_pow3 [r,z]      : axi-symmetric radial velocity deficit ^3 for power estimation [-]
             TI_axial_sym [r,z]           : axi-symmetric radial total turbulence intensity [-]
             TI_meand_axial_sym [r,z]     : axi-symmetric turbulence intensity contribution of wake meandering  [-]
             TI_tot_axial_sym [r,z]       : axi-symmetric total turbulence intensity [-]
             WS [1]                       : rotor averaged wind speed [m/s]
             TI [1]                       : rotor averaged turbulence intensity [-]
             Power [1]                    : available aerodynamic power at turbine from Eq [15] [kW]
             BEM_power [1]                : Blade Element momentum predicted power [kW]
             WS_DWM_BC [1]                : wake averaged wind speed [m/s]
             TI_DWM_BC [1]                : wake averaged turbulence intensity [-]


    """
    dwm = Outputs()
    ###############  Store Flow field #################################################################################
    dwm.WS_axial_ffor          = ffor.WS_axial_ffor
    # dwm.WS_axial_ffor_pow3     = ffor.WS_axial_ffor_pow3
    dwm.TI_axial_ffor          = ffor.TI_axial_ffor
    dwm.TI_meand_axial_ffor    = ffor.TI_meand_axial_ffor
    # dwm.WS_axial_sym           = ffor.WS_axial_sym
    # dwm.WS_axial_sym_pow3      = ffor.WS_axial_sym_pow3
    # dwm.TI_axial_sym           = ffor.TI_axial_sym
    # dwm.TI_meand_axial_sym     = ffor.TI_meand_axial_sym
    # dwm.TI_tot_axial_sym       = np.sqrt(ffor.TI_axial_sym**2 + ffor.TI_meand_axial_sym**2)
    dwm.TI_tot_axial_ffor      = np.sqrt(ffor.TI_axial_ffor**2 + ffor.TI_meand_axial_ffor**2)
    # Store Mean turbine data
    dwm.WS                     = meta.mean_WS_rot
    dwm.TI                     = meta.mean_TI_rot
    dwm.WS_DWM_BC              = meta.mean_WS_DWM
    dwm.TI_DWM_BC              = meta.mean_TI_DWM
    dwm.C2C                    = meta.C2C
    dwm.x_vec                   =  ffor.x_vec
    dwm.y_vec                   =  ffor.y_vec
    dwm.z_vec                   =  ffor.z_vec
    dwm.x_mat                   =  ffor.x_mat
    dwm.y_mat                   =  ffor.y_mat
    dwm.z_mat                   =  ffor.z_mat

    dwm.U_init= mfor.U_init
    dwm.U_init_raw= mfor.U_init_raw
    dwm.vr_mixl = meta.vr_mixl
    dwm.vz_mixl = meta.vz_mixl

    dwm.dA=aero.dA
    dwm.mean_a=aero.mean_a
    # Uniform expansion
    # print aero.mean_a
    # print (1.0-aero.mean_a) / (1.0- ((1.0+meta.fR) * aero.mean_a))
    dwm.f_w=aero.f_w
    dwm.r_w=aero.r_w
    # Boundary conditions for U_w
    dwm.U_w=aero.U_w


    dwm.CP=BEM.CP
    dwm.CPloc=BEM.CPloc
    dwm.CQ=BEM.CQ
    dwm.CQlocCP=BEM.CQloc
    dwm.CPloc=BEM.CPloc
    dwm.CQ=BEM.CQ
    dwm.CQloc=BEM.CQloc
    dwm.CT=BEM.CT
    dwm.CTloc=BEM.CTloc
    dwm.Cd=BEM.Cd
    dwm.Cl=BEM.Cl
    dwm.Cn=BEM.Cn
    dwm.Ct=BEM.Ct
    dwm.Edge=BEM.Edge
    dwm.F=BEM.F
    dwm.Flap=BEM.Flap
    dwm.Fperf=BEM.Fperf
    dwm.Fshen=BEM.Fshen
    dwm.Gamma=BEM.Gamma
    dwm.PITCH=BEM.PITCH
    dwm.PITCH_opt=BEM.PITCH_opt
    dwm.Pn=BEM.Pn
    dwm.Power=BEM.Power
    dwm.Pt=BEM.Pt
    dwm.R=BEM.R
    dwm.RPM=BEM.RPM
    dwm.RPM_opt=BEM.RPM_opt
    dwm.Re=BEM.Re
    dwm.Set=BEM.Set
    dwm.ThrLoc=BEM.ThrLoc
    dwm.ThrLocLn=BEM.ThrLocLn
    dwm.Thrust=BEM.Thrust
    dwm.Torque=BEM.Torque
    dwm.TqLoc=BEM.TqLoc
    dwm.TqLocLn=BEM.TqLocLn
    dwm.Un=BEM.Un
    dwm.Ut=BEM.Ut
    dwm.Vrel=BEM.Vrel
    dwm.a=BEM.a
    dwm.a_last=BEM.a_last
    dwm.alpha=BEM.alpha
    dwm.aprime=BEM.aprime
    dwm.aprime_last=BEM.aprime_last
    dwm.nIt=BEM.nIt
    dwm.phi=BEM.phi
    dwm.r=BEM.r
    dwm.uia=BEM.uia
    dwm.uit=BEM.uit
    dwm.omega=BEM.omega


    dwm.Shear_add_du_dr=mfor.Shear_add_du_dr
    dwm.Shear_add_du_dz=mfor.Shear_add_du_dz
    dwm.TI_DWM=mfor.TI_DWM
    dwm.Turb_Stress_DWM=mfor.Turb_Stress_DWM
    dwm.U=mfor.U
    dwm.U_init=mfor.U_init
    dwm.U_init_pow=mfor.U_init_pow
    dwm.U_init_raw=mfor.U_init_raw
    dwm.V=mfor.V
    dwm.WakeW=mfor.WakeW
    dwm.du_dr_DWM=mfor.du_dr_DWM
    dwm.du_dr_tot=mfor.du_dr_tot
    dwm.visc=mfor.visc

    # print dir(BEM)

    ################### Save to DWM list ###############################################################
    DWM[str(meta.wtg_ind[0])]=dwm

    # if par.get('plot_velocity') is True:
        # if i_wtg==0:
        #     meta.z_vec_old=[]
        # DWM_plot(meta,DWM,i_wtg)

    return DWM


########################################################################################################################
# ***********************************************Dynamic sDWM**********************************************************
# *******************************Working with Mann/LES box approach of meandering**************************************
########################################################################################################################

def get_Meandering_dynamic(meta, meand):


    meand.WakesCentersLocations_in_time = np.load(
            'C:/Users/augus/Documents/Stage/Codes/Mann_Turbulence/Result/Center_Position_in_time_Lillgrund/z_time_center_location.NPY')[meta.iT:]

    # /!\ a mieux placer dans le code /!\
    meand.time = meand.WakesCentersLocations_in_time[0][:, 0]; print 'meand time : ', meand.time
    meand.nt = len(meand.time); print 'number of time points: ', meand.nt
    meta.nt = meand.nt


    return meta, meand


def DWM_MFOR_to_FFOR_dynamic(mfor, meta, meand, ffor):
    """
    Function that calculate the velocity in the fixed (global) frame of reference from the Mfor

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
    plot_bool = True

    ##############################################################################################################
    # recalculate into Cartesian grid
    # initiate/reset Cartesian flow field

    raw_input("Entry to process to MFOR to FFOR")
    DATA_from_Meandering_part = meand.WakesCentersLocations_in_time


    print 'Performing MFoR to FFoR Computation'
    ffor.ffor_flow_field_TI_tmp_tmp = meta.TI * np.ones((meta.nx, meta.ny))  # X = lateral ,Y = vertical
    ffor.TI_axial_ffor_tmp = np.zeros((meta.nx, meta.ny, meta.nz, meta.nt))  # X = lateral ,Y = vertical, time ,Z = streamwise, t = time
    ffor.WS_axial_ffor_tmp = np.zeros((meta.nx, meta.ny, meta.nz, meta.nt))  # X = lateral ,Y = vertical, time ,Z = streamwise, t = time
    ffor.ffor_flow_field_ws_tmp2 = np.zeros((meta.nx, meta.ny, meta.nz, meta.nt))  # X = lateral ,Y = vertical, time ,Z = streamwise, t=time

    ffor.TI_meand_axial_ffor = np.zeros((meta.nx, meta.ny, meta.nz))
    ffor.WS_axial_ffor = np.zeros((meta.nx, meta.ny, meta.nz, meta.nt))
    ffor.TI_axial_ffor = np.zeros((meta.nx, meta.ny, meta.nz, meta.nt))

    ffor.x_vec_t = np.zeros((meta.nx, meta.nz))
    ffor.x_mat_t = np.zeros((meta.nx, meta.ny, meta.nz))
    # Define radial distance vectors
    # r_dist_2              = 999*np.ones((meta.nx,meta.ny))
    r_dist = 999 * np.ones((meta.nx, meta.ny, meta.nt))

    # CREATES THE GLOBAL FLOW FIELDS IN CARTESIAN GRID
    meta.x_mat = np.tile(meta.x_vec.reshape(len(meta.x_vec),1),meta.ny)
    meta.y_mat = np.tile(meta.y_vec,(meta.nx,1))
    meta.z_mat = np.tile(meta.vz, (meta.nx, 1))
    print 'meta vz: ', meta.vz

    # Store the ffor flow field
    ffor.x_vec = (meta.x_vec - meta.hub_x[0]) / 2.
    ffor.y_vec = (meta.y_vec - meta.hub_y) / 2.
    ffor.z_vec = meta.z_vec + np.hstack((0., np.cumsum(meta.hub_z[0:])))[0]
    ffor.x_mat = meta.x_mat / 2.
    ffor.y_mat = meta.y_mat / 2.
    ffor.z_mat = meta.z_mat / meta.dz
    ffor.time = meand.time

    for i_z in np.arange(0, meta.nz, 1):
        # print 'meta.vz[iz]: ', meta.vz[i_z]
        # EXTRACT TI_DWM AND WS_DWM IN MFoR
        # Plot deficit MFOR
        """
        plt.plot(meta.vr_mixl,mfor.U[meta.vz[i_z], :])
        plt.title('Deficit in MFoR at Turbine '+str(7-i_z)+' Location')
        plt.xlabel('vr (polar discretization)')
        plt.ylabel('U (axial velocitie in MFoR)')
        plt.show()
        """

        # After the first Iteration DWM WS DATA need a time dimension, because Ainslie is computed for deficits in time
        #
        try:
            DWM_WS_DATA = mfor.U[meta.vz[i_z], :]

            # print 'DWM_WS_DATA: ', DWM_WS_DATA


        except:
            print 'Fatal error, possibly due to a too low turbulence intensity with respect to the demanded mean wind speed, try increase the input TI'
        DWM_TI_DATA = mfor.TI_DWM[meta.vz[i_z], :]
        # print 'DWM_TI_DATA: ', DWM_TI_DATA

        # Plot TI MFOR

        """
        plt.plot(meta.vr_mixl,mfor.TI_DWM[meta.vz[i_z], :])
        plt.title('TI in MFoR at Turbine '+str(7-i_z)+' Location')
        plt.xlabel('vr (polar discretization)')
        plt.ylabel('TI (Turbulence Intensity in MFoR)')
        plt.show()
        """

        ### Correct DWM_TI_DATA so that no point to have lower TI than "TIamb"
        DWM_TI_DATA[DWM_TI_DATA < np.nanmean(meta.mean_TI_DWM)] = np.nanmean(meta.mean_TI_DWM)

        if plot_bool:
            plt.ion()
            plt.figure(1)
            plt.title('Wake center draw in time')

        for i_t in np.arange(0, len(meand.time), 1):
            #Ro_x = meand.meand_pos_x[i_z, i_t]
            #Ro_y = meand.meand_pos_y[i_z, i_t]
            Ro_x = DATA_from_Meandering_part[i_z][i_t, 1]
            Ro_y = DATA_from_Meandering_part[i_z][i_t, 2]

            #print '(Ro_x,Ro_y): ', [Ro_x, Ro_y]

            #print 'meta.x_mat: ', meta.x_mat
            #r_dist = np.sqrt((meta.x_mat - Ro_x-1.5) ** 2 + (meta.y_mat - Ro_y-1.5) ** 2) # Originally
            r_dist = np.sqrt((ffor.x_mat - Ro_x -0.75) ** 2 + (ffor.y_mat - Ro_y -0.75) ** 2)
            #print 'r_dist: ', r_dist
            # print 'r_dist: ', r_dist

            ###############################################################
            # do we have to keep this for dynamic, or use GLarsen function?
            #print 'mfor.wakeW[meta.vz[i_z]]: ', mfor.WakeW[meta.vz[i_z]]
            tmp_index = r_dist < mfor.WakeW[meta.vz[i_z]] * 1.5
            ################################################################
            #print 'tmp_index: ', tmp_index

            tmp_field_WS = np.ones((meta.nx, meta.ny))
            # print 'tmp_index: ', tmp_index

            # it's here that we change the velocity to be in FFOR
            tmp_field_WS[tmp_index] = np.interp(r_dist[tmp_index], meta.vr_m, DWM_WS_DATA)

            """
            plt.figure()
            plt.subplot(121)
            plt.title('DWM_WS_DATA')
            plt.plot(DWM_WS_DATA, label='DWM_WS_DATA')
            plt.legend()
            plt.subplot(122)
            plt.title('tmp_field_ws')
            plt.plot(tmp_field_WS)
            plt.legend()
            plt.show()
            #"""

            ffor.WS_axial_ffor[:, :, i_z, i_t] = (tmp_field_WS)
            ffor.ffor_flow_field_ws_tmp2[:, :, i_z, i_t] = (tmp_field_WS ** 2)

            tmp_field_TI = meta.TI * np.ones((meta.nx, meta.ny))
            tmp_field_TI[tmp_index] = np.interp(r_dist[tmp_index], meta.vr_m, DWM_TI_DATA)

            ffor.ffor_flow_field_TI_tmp_tmp[:, :] = tmp_field_TI
            ffor.TI_axial_ffor[:, :, i_z, i_t] = ffor.ffor_flow_field_TI_tmp_tmp ** 2



    if plot_bool:
        plt.ion()
        plt.figure()

        x = ffor.x_vec
        y = ffor.y_vec
        X, Y = np.meshgrid(x, y)

        for i_z in np.arange(0, meta.nz):
            for i_t in np.arange(0, meand.nt, 1):
                plt.cla()
                plt.clf()
                plt.subplot(121)
                CS1 = plt.contourf(X, Y, ffor.WS_axial_ffor[:, :, i_z, i_t], 15)
                plt.xlabel('x'), plt.ylabel('y'), plt.title('WS FFoR for Turbine ' + str(i_z)) #7-iz
                plt.colorbar(CS1)

                plt.subplot(122)
                CS1 = plt.contourf(X, Y, ffor.TI_axial_ffor[:, :, i_z, i_t], 15)
                plt.xlabel('x'), plt.ylabel('y'), plt.title('TI FFoR for Turbine ' + str(i_z))  # 7-iz
                plt.colorbar(CS1)
                plt.draw()
                plt.pause(0.001)

        plt.ioff()

    return mfor, ffor, meta, meand


def DWM_get_deficit_FFOR_dynamic(ffor, meta,deficits,ID_waked,inlets_ffor,inlets_ffor_deficits):
    ###############################################################
    DEFI = []
    if meta.steadyBEM_AINSLIE:
        deficits_in_time =  {'1': [], '0': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': []}

        inlets_ffor_in_time =  {'1': [], '0': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': []}
        inlets_ffor_deficits_in_time =  {'1': [], '0': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': []}



    for i_z in np.arange(0, meta.nz, 1):
        deficits_tmp = []
        inlets_ffor_deficits_tmp = []
        inlets_ffor_tmp = []

        # on global frame mesh
        X, Y = np.meshgrid(ffor.x_vec, ffor.y_vec)
        index_trapz = np.sqrt((X + meta.C2C[i_z] / (2. * meta.WTG_spec.R)) ** 2 + (Y) ** 2) >= 0.5

        for i_t in np.arange(0, meta.nt, 1):
            wakedefmask = np.ma.array(np.squeeze(ffor.WS_axial_ffor[:, :, i_z, i_t]), mask=index_trapz,
                                      fill_value=0.0).filled()
            wakedefmasknancoarse = np.ma.array(np.squeeze(ffor.WS_axial_ffor[:, :, i_z, i_t]), mask=index_trapz,
                                               fill_value=np.nan).filled()
            disk = np.ma.array(np.zeros(wakedefmask.shape), mask=~index_trapz, fill_value=1.0).filled()
            disk_area = np.trapz(np.trapz(disk, dx=1. / meta.dy), dx=1. / meta.dx)
            trapz2 = np.trapz(np.trapz(wakedefmask, dx=1. / meta.dy), dx=1. / meta.dx)

            # on finer mesh
            # values=np.squeeze(ffor.WS_axial_ffor[:,:,i_z])
            # X, Y = np.meshgrid(ffor.x_vec,ffor.y_vec)
            # points=np.vstack((np.ravel(X),np.ravel(Y)))
            # grid_x, grid_y = np.mgrid[-0.5:0.5:meta.dR*36j, -0.5:0.5:meta.dR*36j]
            # wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='cubic')
            #
            # X2,Y2=np.meshgrid(np.linspace(-0.5,0.5,36*meta.dR),np.linspace(-0.5,0.5,36*meta.dR))
            # index_trapz=np.sqrt((X2 + meta.C2C[i_z]/(2.*meta.WTG_spec.R))**2 + (Y2)**2 )>=0.5
            #
            # wakedefmask = np.ma.array(wake_i, mask=index_trapz, fill_value=0.0).filled()
            # # wakedefmasknan = np.ma.array(wake_i, mask=index_trapz, fill_value=np.nan).filled()
            # disk = np.ma.array(np.zeros(wakedefmask.shape), mask=~index_trapz, fill_value=1.0).filled()
            # disk_area=simps(simps(disk,np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dy)),np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dx))
            # trapz2=simps(simps(wakedefmask,np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dy)),np.linspace(-0.5,0.5,36*meta.dR),dx=1./(36.*meta.dx))
            #
            # COMPUTING TO PLOT
            DEFI = DEFI + [trapz2 / disk_area]

            deficits_tmp.append(trapz2 / disk_area)
            inlets_ffor_deficits_tmp.append(wakedefmasknancoarse)
            inlets_ffor_tmp.append([np.vstack(
                ((meta.x_vec - meta.hub_x[i_z]), (meta.y_vec - meta.hub_y), ffor.WS_axial_ffor[:, :, i_z, i_t]))])
            #print 'deficits_tmp: ', trapz2 / disk_area

        if meta.steadyBEM_AINSLIE:  # see meta definition in cDWM
            # think a way to save the dynamic DATA
            deficits_in_time[str(meta.wtg_ind[i_z])].append(deficits_tmp)
            inlets_ffor_deficits_in_time[str(meta.wtg_ind[i_z])].append(inlets_ffor_deficits_tmp)
            inlets_ffor_in_time[str(meta.wtg_ind[i_z])].append(inlets_ffor_tmp)
            print 'inlets_ffor_tmp shape: ', np.shape(inlets_ffor_tmp)
            print 'inlets_ffor_deficits_tmp shape: ', np.shape(inlets_ffor_deficits_tmp)
            # Average in time
            deficits_tmp = np.mean(deficits_tmp)
            inlets_ffor_deficits_tmp = np.mean(inlets_ffor_deficits_tmp, axis=0)
            inlets_ffor_tmp = np.mean(inlets_ffor_tmp, axis=0)
            print 'inlets_ffor_deficits shape (after mean): ', np.shape(inlets_ffor_tmp)
            print 'inlets_ffor_deficits_tmp shape (after mean): ', np.shape(inlets_ffor_deficits_tmp)

        # Updating

        deficits[str(meta.wtg_ind[i_z])].append(deficits_tmp)
        inlets_ffor_deficits[str(meta.wtg_ind[i_z])].append(inlets_ffor_deficits_tmp) #conatained NaN value
        inlets_ffor[str(meta.wtg_ind[i_z])].append(inlets_ffor_tmp)

        ID_waked[str(meta.wtg_ind[i_z])].append(meta.wtg_ind[0])
    print 'deficits dict: ', deficits
    """
    if not meta.steadyBEM_AINSLIE:
        i_z = 1
        plt.figure()
        plt.title('WT' + str(i_z) + ' deficits in time')
        plt.plot(range(0, meta.nt),deficits[str(meta.wtg_ind[i_z])][0])
        plt.xlabel('Simulation time iteration, Nt'), plt.ylabel('Deficit')
        plt.show()
    if meta.steadyBEM_AINSLIE:
        i_z = 1
        plt.figure()
        plt.title('WT' + str(i_z) + ' deficits in time')
        plt.plot(range(0, meta.nt), deficits_in_time[str(meta.wtg_ind[i_z])][0], label='temporal data')
        plt.plot(range(0,meta.nt), [deficits[str(meta.wtg_ind[i_z])][0] for i_t in range(0, meta.nt)], label='Average deficit in time')
        plt.xlabel('Simulation time iteration, Nt'), plt.ylabel('Deficit'), plt.legend()
        plt.show()
    #"""
    # Plotting
    """

        print trapz2 / disk_area
    plt.plot(np.arange(0,meta.nz,1),DEFI)
    plt.title('Average Deficits for each Turbine'),plt.xlabel('Turbine'),plt.ylabel('Deficits')
    plt.show()
        #####Print
        #print deficits
        #print inlets_ffor_deficits
        #print inlets_ffor
    """
    """
        print wakedefmasknancoarse
    plt.plot(np.arange(0,meta.nz,1),DEFI)
    plt.title('Deficits for each Turbine'),plt.xlabel('Turbine'),plt.ylabel('Deficits')
    plt.show()
        #####Print
        #print deficits
        #print inlets_ffor_deficits
        #print inlets_ffor
    """
    raw_input('Press any key to continue')
    return deficits, ID_waked, inlets_ffor, inlets_ffor_deficits


def DWM_get_turb_dynamic(ffor,meta,turb,inlets_ffor_turb,):
    """
    Function that calculate the rotor turbulence intensity in time
    Parameters
    ----------
    ffor:(instance of class): Instance of class Ffor holding the fixed frame of reference velocity field in global
    WF coordinates
    meta (instance of class): Instance of class Meta holding DWM core variables
    turb:  dict(nWT) holding a list of turbulence intensities contributions from upstream wakes
    inlets_ffor_turb: dict(nWT) holding a list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position

    Returns
    -------
    turb:  dict(nWT) updated list of turbulence intensities contributions from upstream wakes
    inlets_ffor: dict(nWT) updated list of array containing the flow field in the fixed frame of reference from upstream wakes contributions
    inlets_ffor_turb: dict(nWT) updated list of array containing the turbulence field in the fixed frame of reference from upstream wakes contributions at the rotor position
    """
    PLOT=[]
    if meta.steadyBEM_AINSLIE:
        turb_in_time =  {'1': [], '0': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': []}
        inlets_ffor_turb_in_time =  {'1': [], '0': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': []}
    for i_z in np.arange(0,meta.nz,1):
        X, Y = np.meshgrid(ffor.x_vec, ffor.y_vec)
        index_trapz = np.sqrt((X + meta.C2C[i_z] / (2. * meta.WTG_spec.R)) ** 2 + (Y) ** 2) >= 0.5

        turb_tmp = []
        inlets_ffor_turb_tmp =[]
        for i_t in np.arange(0, meta.nt, 1):

            turbmask = np.ma.array(np.squeeze(ffor.TI_axial_ffor[:, :, i_z, i_t]), mask=index_trapz, fill_value=0.0).filled()
            turbmasknan = np.ma.array(np.squeeze(ffor.TI_axial_ffor[:, :, i_z, i_t]), mask=index_trapz, fill_value=np.nan).filled()
            disk = np.ma.array(np.zeros(turbmask.shape), mask=~index_trapz, fill_value=1.0).filled()
            disk_area=np.trapz(np.trapz(disk,dx=1./meta.dy),dx=1./meta.dx)
            trapz2=np.trapz(np.trapz(turbmask,dx=1./meta.dy),dx=1./meta.dx)

            turb_tmp.append(trapz2/disk_area)
            inlets_ffor_turb_tmp.append(turbmasknan)

            #print 'trapz2/disk_area: ', trapz2/disk_area

        if meta.steadyBEM_AINSLIE:   # see meta definition in cDWM
            # think a way to save the dynamic DATA
            turb_in_time[str(meta.wtg_ind[i_z])].append(turb_tmp)
            inlets_ffor_turb_in_time[str(meta.wtg_ind[i_z])].append(inlets_ffor_turb_tmp)
            print 'inlets_ffor_turb_tmp shape: ', np.shape(inlets_ffor_turb_tmp)
            # Average in time
            turb_tmp = np.mean(turb_tmp)
            inlets_ffor_turb_tmp = np.mean(inlets_ffor_turb_tmp, axis=0)
            print 'inlets_ffor_turb_tmp shape (after mean): ', np.shape(inlets_ffor_turb_tmp)
        turb[str(meta.wtg_ind[i_z])].append(turb_tmp)
        inlets_ffor_turb[str(meta.wtg_ind[i_z])].append(inlets_ffor_turb_tmp)

    #PLOTTING
    """
    plt.plot(np.arange(0,meta.nz,1),PLOT)
    plt.title('Average Turbulence for each Turbine'),plt.xlabel('Turbine Location'), plt.ylabel('TI')
    plt.show()
    #"""
    """
    i_z = 1
    plt.figure()
    plt.title('WT' + str(i_z) + ' turb in time')
    plt.plot(range(0, meta.nt), turb[str(meta.wtg_ind[i_z])][0])
    plt.show()
    #"""
    raw_input('End of get turb process (press any key to continue)')
    return turb,inlets_ffor_turb


def DWM_aero_dynamic(meta,ffor,aero,deficits,turb,inlets_ffor,inlets_ffor_deficits,out,ID_waked):
    """ Aerodynamique module of the DWM. This module contains the wake summation module (deficit and turbulence accumulation)
    The steady state blade element momentum

        Inputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        aero (instance of class): Instance of class Aero holding BEM-aero core variables

        deficits: dict(nWT) holding a list of deficits contributions from upstream wakes
        turb:  dict(nWT) holding a list of mean turbulence intensities contributions from upstream wakes
        inlets_ffor: dict(nWT) holding a list of array containing the flow field in the fixed frame of reference from
        upstream wakes contributions
        inlets_ffor_deficits: dict(nWT) holding a list of array containing the flow field in the fixed frame of
        reference from upstream wakes contributions at the rotor position
        out: dict(nWT),holding the main sDWM model outputs i.e mean power from BEM, mean power estimated from powercurve,
        mean rotor averaged wind speed, mean rotor average turbulence intensity, mean thrust coefficient from BEM and
        from power curve

        Outputs
        ----------
        aero (instance of class): updated Instance of class Aero holding BEM-aero core variables
        mfor (instance of class): updated Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
        out dict(nWT): dict including mean power from PC and BEM, mean thrust coefficient from PC and BEM
        BEM (instance of class): holds the key results from the BEM calculation

    """
    mfor   =  MFoR(meta.WTG)
    ## Compute the average wake deficit accumulation

    if meta.wake_ind_setting == 1:
        meta.mean_WS_DWM = []
        Deficits_in_time = deficits.get(str(meta.wtg_ind[0]))  # we get the deficit of the generating plan for the current WT
        print 'Deficits_in_time for the current WT: ', Deficits_in_time
        if not Deficits_in_time: # means that the current turbine is in free stream (first Turbine in the row)
            nt=1  # for the first Turbine in the free stream, there is no temporal iteration

        else:
            nt=meta.nt
            Deficits_in_time = np.array(Deficits_in_time)
        print 'Deficits_in_time for the current WT: ', Deficits_in_time
        print 'Deficits_in_time for the current WT shape: ', np.shape(Deficits_in_time)
        raw_input('Press Enter to continue')
        for i_t in np.arange(0, nt, 1):
            if nt == 1:   # means that the current turbine is in free stream (first Turbine in the row)
                rWS = np.array([1.0])
            else :
                rWS = Deficits_in_time[:, i_t]
            print 'deficit of the generating plan for the current WT at a certain time: ', rWS

            if meta.accu == 'linear':
                meta.mean_WS_DWM.append(meta.WS*(1.-(np.sum([1. - xx for xx in rWS]))))
            elif meta.accu == 'quadratic':
                meta.mean_WS_DWM.append(meta.WS*(1.-np.sqrt(np.sum([(1.-xx)**2 for xx in rWS]))))
            elif meta.accu == 'dominant':
                meta.mean_WS_DWM.append(meta.WS*(1.-(np.max([1. - xx for xx in rWS]))))
                # meta.mean_WS_DWM= meta.WS*(1.-(np.max([1. - xx for xx in rWS])))
            # elif meta.accu == 'bypassed':
            #     meta.mean_WS_DWM= meta.WS
            elif meta.accu == 'ewma':
                print 'This model is currently being developped and implemented soon'
            else:
                print 'You have not specified any wake accumulation procedure... or mispelled it. Options are ''linear'', ''quadratic'' , ''dominant'' or ''ewma'''
    else:
        meta.mean_WS_DWM   = meta.WS # free stream everywhere similar to current implementation of HAWC2
    meta.mean_WS_DWM = np.array(meta.mean_WS_DWM)
    print 'meta.mean_WS_DWM: ', meta.mean_WS_DWM
    # Set buildup of turbulence
    if meta.Tbuildup_setting==1:
        ti = turb.get(str(meta.wtg_ind[0]))

        Turb_in_time = deficits.get(str(meta.wtg_ind[0]))  # we get the deficit of the generating plan for the current WT
        print 'Turb_in_time for the current WT: ', Turb_in_time
        if not Turb_in_time:  # means that the current turbine is in free stream (first Turbine in the row)
            nt = 1  # for the first Turbine in the free stream, there is no temporal iteration

        else:
            nt = meta.nt
            Turb_in_time = np.array(Turb_in_time)
        print 'Turb_in_time for the current WT: ', Turb_in_time
        print 'Turb_in_time for the current WT shape: ', np.shape(Turb_in_time)
        raw_input('Press Enter to continue')
        for i_t in np.arange(0, nt, 1):

            if nt == 1:
                meta.mean_TI_DWM  = meta.TI
            else:
                meta.mean_TI_DWM  = np.max(Turb_in_time[:, i_t])
    else:
        meta.mean_TI_DWM  = meta.TI

    meta.mean_TI_DWM = np.array(meta.mean_WS_DWM)
    print 'meta.mean_TI_DWM: ', meta.mean_TI_DWM

    raw_input('Run BEM at accumulated deficit')
    # Run BEM at accumulated deficit
    aero,BEM   =  DWM_rotor_aero(meta,aero,ID_waked)

    # domain induction
    a_domain     = np.interp(meta.vr_m,np.hstack(([aero.r_w, aero.r_w[-1]+0.01, aero.r_w[-1]+0.02])), np.hstack((( 1.0 - aero.U_w), [0., 0.])))

    ## Compute the accumulated flow field for accurate inlet definition of the MFoR wake calculation
    # center all disks before wake accumulation
    if not inlets_ffor.get(str(meta.wtg_ind[0])):# if turbine in FREESTREAM, then mfor.Uinit is initialized in the calc_mixL module
        # mfor.U_init =   None # set to none for further initialization
        radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        mfor.U_init=f3(meta.vr_mixl)
        mfor.U_init=smooth( mfor.U_init,window_len=5)
    elif meta.accu_inlet is False:
        # mfor.U_init =   None
        radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        mfor.U_init=f3(meta.vr_mixl)
        mfor.U_init=smooth( mfor.U_init,window_len=5)
    else:   # if turbine not in the freestream, we need to compute the proper accumulated inlet to the turbine
        ranger=np.linspace(-1.,1.,meta.dR*2.)  # np.linspace(-2.,2.,meta.dR*4.)
        inlets_ffor_deficits_np_3D=np.ones((len(ranger) ,len(ranger)   , len(inlets_ffor[str(meta.wtg_ind[0])])))
        # grid_x, grid_y = np.mgrid[-1.:1.:meta.dR*2j, -1.:1.:meta.dR*2j]
        grid_x, grid_y = np.mgrid[-1.:1.:meta.dR*2j, -1.:1.:meta.dR*2j]
        for ii in range(len(inlets_ffor[str(meta.wtg_ind[0])])):
            offsets=2.*(min(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0])+abs(max(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0])-min(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0]))/2.)
            # need to interp on a new array of equal size
            values=inlets_ffor_deficits[str(meta.wtg_ind[0])][ii]
            X, Y = np.meshgrid(inlets_ffor[str(meta.wtg_ind[0])][ii][0][0]-offsets,inlets_ffor[str(meta.wtg_ind[0])][ii][0][1])
            points=np.vstack((np.ravel(X),np.ravel(Y)))
            # wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='linear')
            wake_i=interpolate.griddata(points.T,np.ravel(values),(grid_x, grid_y), method='linear')
            inlets_ffor_deficits_np_3D[:,:,ii]=wake_i
            # print wake_i
        if meta.accu == 'linear':
            U_init=1.-(np.sum(1.-inlets_ffor_deficits_np_3D,axis=2))
        elif meta.accu == 'quadratic':
            U_init=1.-np.sqrt(np.sum((1.-inlets_ffor_deficits_np_3D)**2,axis=2))
        elif meta.accu == 'dominant':
            U_init=np.amin(inlets_ffor_deficits_np_3D, axis=2)
        elif meta.accu == 'ewma':
            print 'This model is currently being developped and implemented soon'
        else:
            print 'You have not specified any wake accumulation procedure... or mispelled it. Options are ''linear'', ''quadratic'' , ''dominant'' or ''ewma'''
        # Transform to axisymmetric profile inlet
        r_dist_2= np.sqrt(grid_x**2 + grid_y**2 )  #Find distance to centre of wake plane
        ffor.WS_axial_sym      = np.ones((len(np.arange(0,meta.dR+1.,1))))
        ffor.WS_axial_sym[0]=np.nanmean(U_init[r_dist_2 < (1.05*np.amin(r_dist_2))])
        print 'meta.dR: ', meta.dR
        for i_r_pos in np.arange(1,meta.dR+1,1):
            print 'i_r_pos: ', i_r_pos
            a=r_dist_2 > ((i_r_pos+1-1.5)*(1.0/meta.dR))# rotor neg boundaries
            bb=r_dist_2 < ((i_r_pos+1-0.5)*(1.0/meta.dR)) #rotor pos boundaries
            c=np.logical_and(a,bb)
            bin_filter = c
            tmp_ffor_flow_field_ws_mean          = U_init[bin_filter]
            ffor.WS_axial_sym[i_r_pos]           = np.nanmean(tmp_ffor_flow_field_ws_mean)
        ffor.r_sym=np.arange(0,meta.dR+1.,1)/meta.dR
        # Update the DWM inlet
        if ffor.r_sym[-1] >= meta.vr_mixl[-1]:
            # mfor.U_init = (1.0-a_domain) * np.interp(meta.vr_mixl,ffor.r_sym,ffor.WS_axial_sym)
            print 'spline interpolation'
            mfor.U_init = (1.0-a_domain) * interpolate.InterpolatedUnivariateSpline(meta.vr_mixl,ffor.r_sym,ffor.WS_axial_sym)
        else:
            print 'meta.lr_mixl: ', meta.lr_mixl
            mfor.U_init = (1.0-a_domain) * np.hstack((ffor.WS_axial_sym.ravel(), np.ones((((meta.dR * int(meta.lr_mixl))-ffor.WS_axial_sym.size),1)).ravel()))
        # Finishing
        mfor.U_init_raw=mfor.U_init
        mfor.U_init[mfor.U_init < 0.0]=0.0 # prevent from negative velocities on linear summation
        mfor.U_init=smooth( mfor.U_init,window_len=5)

    # Power curve based
    try:
        if BEM.derated is False:
            # print 'use standard ws for curve'
            aero.pow_cur=meta.WTG_spec.get_P(meta.mean_WS_DWM)
            aero.ct_cur=meta.WTG_spec.get_CT(meta.mean_WS_DWM)
        else:
            # print 'use demanded ws for curve'
            aero.pow_cur=meta.WTG_spec.get_P(BEM.Ud)
            aero.ct_cur=meta.WTG_spec.get_CT(BEM.Ud)
    except:
        aero.pow_cur=0.
        aero.ct_cur=0.
    # write outlets

    #/!\/!\ not put in commentary this  /!\/!\
    """
    out[str(meta.wtg_ind[0])]=[]
    out[str(meta.wtg_ind[0])].append(float(format(aero.Power/1000., '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(meta.mean_WS_DWM, '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(meta.mean_TI_DWM, '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(aero.CT/1., '.2f')))
    out[str(meta.wtg_ind[0])].append(float(format(aero.pow_cur, '.2f'))) # based on power curve
    out[str(meta.wtg_ind[0])].append(float(format(aero.ct_cur, '.2f'))) # based on power curve
    #"""
    return aero, mfor, out, BEM