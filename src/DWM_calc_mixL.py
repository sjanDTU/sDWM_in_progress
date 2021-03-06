# -*- coding: utf-8 -*-
""" Python interpretation of Rolf-Erik Keck's Ainslie model
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""

import numpy as np
from math import pi
from scipy import io, interpolate, linalg
from DWM_misc import smooth
import matplotlib.pyplot as plt


def DWM_init_calc_mixl(meta,aero,mfor):
    """
        Inputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        mfor (instance of class): Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model

        Outputs
        ----------
        mfor (instance of class): Update instance of class Mfor with initialized flow variables based on wake case
        F1_vector (np.array float): DWM filter functions F 1 governing the development of turbulent stresses ref[6-7], eq 3 of [5]
        F2_vector (np.array float): DWM filter functions F 2 governing the development of turbulent stresses ref[6-7], eq 3 of [5]
        visc_wake1 (np.array float): initialized contribution due to ambient turbulence for the eddy viscosity ref[5], Eq 3 first term of right and side
        visc_wake2 (np.array float): initialized contribution due to shear layer of wake deficit for the eddy viscosity ref[5], Eq 3 last term of right and side
        visc_wake (np.array float): eddy viscosity of DWM, combination of visc_wake1 and visc_wake2
        u_star_DEF (np.array float): non dimensionalivelocity scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
        l_star_DEF (np.array float): non dimensional integral length scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
        One_div_du_dr_DWW (np.array float): denominator of Eq (5) in ref [5]
        width (np.array int): wake width vector index
    """

    F1_vector, F2_vector = None , None
     # Sets the input vector
    # bin_filter = np.array([aero.r_w.max(0) > meta.vr_mixl]).astype(int).flatten()
    # xq=meta.vr_mixl[:sum(bin_filter)]*bin_filter[:sum(bin_filter)]
     # new approach
    # bin_filter = np.array([aero.r_w.max(0) > meta.vr_mixl]).astype(bool).flatten()
    # xq=np.hstack((meta.vr_mixl[bin_filter],aero.r_w.max(0)))
    # if mfor.U_init is None:
        # mfor.U_init_raw= np.hstack((np.interp(xq,aero.r_w, aero.U_w),np.linspace(meta.U0,meta.U0,(len(meta.vr_mixl[sum(bin_filter):])))))
        # mfor.U_init = np.hstack((np.interp(xq,aero.r_w, aero.U_w),np.linspace(meta.U0,meta.U0,(len(meta.vr_mixl[sum(bin_filter):])))))

        # new approach
        # mfor.U_initi = np.interp(xq,aero.r_w, aero.U_w)
        # mfor.U_initi=np.hstack((  mfor.U_initi, np.linspace(meta.U0, meta.U0, (len(meta.vr_mixl[sum(bin_filter):])))))
        # vr_mixli=np.hstack(( xq  , meta.vr_mixl[np.invert(bin_filter)] ))
        # mfor.U_init=np.interp(meta.vr_mixl,vr_mixli,mfor.U_initi)
        # mfor.U_initxbh=np.interp(np.linspace(vr_mixli[0],vr_mixli[-1],400),vr_mixli,mfor.U_initi)
        # mfor.U_init=np.interp(meta.vr_mixl,np.linspace(vr_mixli[0],vr_mixli[-1],400), mfor.U_initxbh)

        # Much better approach
        # radial=np.hstack((aero.r_w,aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        # vel=np.hstack((aero.U_w,1., 1.))
        # radial=np.hstack((aero.r_w[0:-2],aero.r_w.max(0)+0.01,meta.vr_mixl[-1]))
        # vel=np.hstack((aero.U_w[0:-2],0.99, 1.))
        # f3=interpolate.InterpolatedUnivariateSpline(radial,vel,  k=1)
        # mfor.U_init=f3(meta.vr_mixl)
        # mfor.U_init=smooth( mfor.U_init,window_len=5)
    # Generate the DWM filter functions for the eddy viscosity formulation

    # ----------------------- # Filter Function Definitions # -------------------------------------------------------- #
    if not meta.without_filter_functions:
        if meta.Keck:
            F1_vector  = np.hstack((np.linspace(meta.f1[0],1,meta.f1[1]*meta.dz/2),np.linspace(1,1,meta.lz_mixl*meta.dz/2)))
            F2_z_vec   = np.arange(2+1./meta.dz,(len(F1_vector)+1)*(1./meta.dz),1./meta.dz)
            F2_vector  = np.hstack((np.linspace(meta.f2[0],meta.f2[0],2*meta.dz), 1.-(1.-meta.f2[0])*np.exp(-meta.f2[1]*(F2_z_vec-2.))))

        if meta.Madsen:
            # Same F1 as Keck
            # Not the same F2 as Keck... I don't find the formula so I interpolate the graph.
            F1_vector = np.hstack((np.linspace(meta.f1[0], 1, meta.f1[1] * meta.dz / 2), np.linspace(1, 1, meta.lz_mixl * meta.dz / 2)))
            F2_vector = meta.F2(meta.vz_mixl / 2)
            # F2 vector seems to be constituted by 4 functions
            # [0D]0.07 -> 0.07 [2D],
            # linear curve 0.07[2D] -> 0.3 [7D],
            # exponential curve  0.3[7D] -> 1[10D]
        if meta.Larsen:
            # We have to introduce the non linear function F_amb
            # I don't found the formula for f_amb I don't find the formula so I interpolate the graph.
            # it seems to be same F2 as Madsen
            k = 2.5
            F1_vector = (np.arctan(k * (meta.vz_mixl - 5)) / pi + 0.5)
            F2_vector = meta.F2(meta.vz_mixl/2)
            meta.F_amb = 0.12/meta.mean_TI_DWM

            # F_amb curve in 1/TI_amb, analytical solution possible => 0.12/TI_amb
            # F1 curve in arctan(x/D)

    if meta.AINSLIE_EV_details:
        print "Initial lenghts of F1 and F2 don't correspond to vz_mixl"
        print "Plot: We restricted the plot to vz_mixl domain"
        plt.figure('Filter Function')
        plt.title('Filter Function for Ainslie Calculations')
        plt.plot(meta.vz_mixl, F1_vector[:len(meta.vz_mixl)], label='F1')
        plt.plot(meta.vz_mixl, F2_vector[:len(meta.vz_mixl)], label='F2')
        plt.xlabel('Dowstream position z (Mixl Domain) [R]'), plt.ylabel('[-]')
        plt.legend(), plt.show()

    # initiate the U, V, visc etc... matrices
    mfor.V                 = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float)
    mfor.U                 = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.visc              = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.du_dr_DWM         = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.du_dr_tot         = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.Turb_Stress_DWM   = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    mfor.TI_DWM            = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    visc_wake         = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    visc_wake1        = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    visc_wake2        = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )
    One_div_du_dr_DWM = np.zeros((len(meta.vz_mixl), len(meta.vr_mixl)),dtype=float )

    #print 'shape(mfor.V)=shape(mfor.U): ', np.shape(mfor.V)      # Example: (306L, 80L)
    # DWM boundary conditions

    # Rotor plane
    mfor.U[0,:] = mfor.U_init
    #print mfor.U[0,:]
    #print 'len(U) in mixl domain: '+str(len(mfor.U[0,:]))

    #Centerline
    mfor.V[0,:]     = 0.0
    #print 'len(V) in mixl domain: ' + str(len(mfor.V[0, :]))
    #print 'len(vr_mixl) polar coordinate: '+str(len(meta.vr_mixl))
    #plt.plot(meta.vr_mixl,mfor.U[0,:])
    #plt.show()

    # Atmospheric stability effects based on Keck et al. [4]
    if (meta.atmo_stab== 'VU')==1:
          L_ABL_vector         = [42.9157,      68.5912,      88.0709]
          UW_UU_vector         = [-0.27991,    -0.26012,    -0.23296]
          L_DEF_vector         = [14.8641,      19.5041,      22.0656]
          UU_DEF_UU_ABL_vector = [0.51678,     0.47861,     0.42581]
          UW_DEF_UU_DEF_vector = [-0.23688,    -0.13927,   -0.097113]
    elif (meta.atmo_stab== 'U')==1:
          L_ABL_vector         = [39.0648,      61.8843,      80.8478]
          UW_UU_vector         = [-0.27441,    -0.28051,    -0.26325]
          L_DEF_vector         = [14.0137,      18.1673,      20.8881]
          UU_DEF_UU_ABL_vector = [0.51067,      0.4433,     0.41874]
          UW_DEF_UU_DEF_vector = [-0.2544,    -0.18689,    -0.12968]
    elif (meta.atmo_stab== 'NU')==1:
          L_ABL_vector         = [33.2038,       48.321,      66.7102]
          UW_UU_vector         = [-0.27136,    -0.28125,    -0.28078]
          L_DEF_vector         = [12.7207,      15.9199,      18.7976]
          UU_DEF_UU_ABL_vector = [0.54998,     0.49571,      0.4139]
          UW_DEF_UU_DEF_vector = [-0.26641,    -0.22075,    -0.18199]
    elif (meta.atmo_stab== 'N')==1:
          L_ABL_vector         = [26.5352,       34.026,      40.7458]
          UW_UU_vector         = [-0.27359,    -0.27887,    -0.27935]
          L_DEF_vector         = [11.065,      12.9746,      14.4395]
          UU_DEF_UU_ABL_vector = [0.63044,     0.57982,      0.5287]
          UW_DEF_UU_DEF_vector = [-0.27341,    -0.25684,    -0.24217]
    elif (meta.atmo_stab== 'NS')==1:
          L_ABL_vector         = [21.2064,      27.1416,      34.2689]
          UW_UU_vector         = [-0.27331,    -0.27636,    -0.28028]
          L_DEF_vector         = [9.50836,      11.2453,      13.0561]
          UU_DEF_UU_ABL_vector = [0.69202,     0.63823,     0.59067]
          UW_DEF_UU_DEF_vector = [-0.27954,    -0.26985,    -0.25258]
    elif (meta.atmo_stab== 'S')==1:
          L_ABL_vector         = [12.4648,      14.6943,      21.5762]
          UW_UU_vector         = [-0.27085,    -0.27464,     -0.2781]
          L_DEF_vector         = [6.3775,      7.2553,      9.6425]
          UU_DEF_UU_ABL_vector = [0.80217,     0.78088,     0.71086]
          UW_DEF_UU_DEF_vector = [-0.28288,     -0.2835,     -0.2758]
    elif (meta.atmo_stab== 'VS')==1:
          L_ABL_vector         = [7.4849,      10.4295,      23.5966]
          UW_UU_vector         = [-0.26527,     -0.2729,    -0.26608]
          L_DEF_vector         = [4.21111,      5.53777,      10.2117]
          UU_DEF_UU_ABL_vector = [0.87381,     0.83775,     0.63848]
          UW_DEF_UU_DEF_vector = [-0.27585,    -0.28304,    -0.27871]
    L_ABL          = np.interp(meta.WTG_spec.H,[40.,  100.,  160.], L_ABL_vector  )    # int_Lww(k3) i.e (Integral length scale (in vertical directions), from ww(k3))
    UW_UU          = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], UW_UU_vector )    # ratio of UW and UU stresses for whole spectra
    L_DEF          = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], L_DEF_vector )   #    Integral length scale (in vertical directions), Meandering length scale subtracted, from ww(k3)
    UU_DEF_UU_ABL  = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], UU_DEF_UU_ABL_vector)    # Part of normal stress in the deficit module
    UW_DEF_UU_DEF  = np.interp( meta.WTG_spec.H,[40.,  100.,  160.], UW_DEF_UU_DEF_vector)    # ratio of UW and UU stres
    # DO NOT CHANGE!!
    Rotor_R        = 40. #ATMOSTAB ANALYSIS IS CARRIED OUT OVER R = 40m, which should be used to normalize the length scales
    l_star_ABL     = L_ABL / Rotor_R#
    l_star_DEF     = L_DEF / Rotor_R#

    # Normalize UU_160m to neutral condition (to use calibration from Keck et al. [3])
    UU_DEF_UU_ABL_fac  = np.interp(meta.WTG_spec.H,[40., 100., 160.], [0.63044,     0.57982,      0.5287])
    UU_DEF_UU_ABL      = UU_DEF_UU_ABL / UU_DEF_UU_ABL_fac

    #CALCULATE u* according to:
    #1. u* ~= (mean(u'w')^2 )^0.25
    #2. {mean(u'w') = mean(u'u')*Cuw_uu}
    #3. {u' ~= TI (in normalized form)}
    # => u* ~= ((TI^2 * Cuw_uu )^2)^0.25
    # u_star_ABL     = ( (  (meta.mean_TI_DWM[0])**2            * np.abs(UW_UU)         )**2 )**0.25 # later on we divide by 100 and not 1000 and TI=10 not 100!!!
    u_star_ABL     = ( (  (meta.TI)**2            * np.abs(UW_UU)         )**2 )**0.25 # replaced by inflow TI
    u_star_DEF     = ( (  (meta.mean_TI_DWM)**2 * UU_DEF_UU_ABL * np.abs(UW_DEF_UU_DEF) )**2 )**0.25
    mfor.Shear_add_du_dz = u_star_ABL / l_star_ABL
    width=np.zeros((len(np.arange(1,len(meta.vz_mixl),1)), 1),dtype=int)
    return meta,mfor,F1_vector,F2_vector,visc_wake1,visc_wake2,visc_wake,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,width


def DWM_calc_wake_width(mfor,width,b_loop,meta,j):
    """Function that estimate the wake width

        Inputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        mfor (instance of class): Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
        width (np.array int) : wake width vector index
        j (int) : index in the main Ainslie forward scheme loop

        Outputs
        ----------
        width (np.array int): updated wake width vector index
    """


    # Calculating wake width
    #print 'meta.dr_mixl: ', meta.dr_mixl
    #print 'meta.lr_mixl: ', meta.lr_mixl
    dr_DWM        = 1.0/meta.dr_mixl
    #print 'dr_DWM = 1.0/meta.dr_mixl: ', dr_DWM
    # meta.lr_mixl=18.95
    r_vec_DWM     = np.arange(dr_DWM/2.0,meta.lr_mixl-dr_DWM/2.0+dr_DWM,dr_DWM)
    #print 'r_vec_DWM = np.arange(dr_DWM/2.0,meta.lr_mixl-dr_DWM/2.0+dr_DWM,dr_DWM): ', r_vec_DWM

    r_vec_DWM = r_vec_DWM[r_vec_DWM<meta.lr_mixl] # bug fix
    #print 'r_vec_DWM = r_vec_DWM[r_vec_DWM<meta.lr_mixl] # bug fix: ', r_vec_DWM

    dA_DWM        = np.array([pi*r_vec_DWM[1:]**2-pi*r_vec_DWM[0:-1]**2])
    #print 'dA_DWM = np.array([pi*r_vec_DWM[1:]**2-pi*r_vec_DWM[0:-1]**2]): ', dA_DWM


    # print meta.lr_mixl
    # print 'r_vec_DWM',r_vec_DWM
    # print 'dA_DWM.shape', dA_DWM.shape
    # # print meta.lr_mixl-dr_DWM/2.0+dr_DWM
    # # print np.linspace(0,meta.lr_mixl-meta.dr,(meta.lr_mixl)/meta.R_WTG*meta.dr_mixl)
    # # print meta.lx
    # print 'len(meta.vr_mixl)', len(meta.vr_mixl)
    # print 'mfor.U.shape', mfor.U.shape
    # print 'np.linspace(0,meta.lr_mixl-meta.dr,(meta.lr_mixl)/meta.R_WTG*meta.dr_mixl).shape', np.linspace(0,meta.lr_mixl-meta.dr,(meta.lr_mixl)/meta.R_WTG*meta.dr_mixl).shape

    Def_DWM       = np.sum(( 1.0-mfor.U[j-1,1:])*dA_DWM  )
    #print 'Def_DWM = np.sum(( 1.0-mfor.U[j-1,1:])*dA_DWM ):'  #compute the deficit containing in each dA_DWM
    #print Def_DWM
    #Def_DWM       = sum((1 - U(j-1,2:end) ).* dA_DWM);
    # if j==1:
    #     tsolve = time.time()
    while(b_loop):
        b_loop = b_loop +1
        Def_DWM_mixL = np.sum((1.0 - mfor.U[j-1,1:(b_loop)]) * dA_DWM[0,0:(b_loop-1)])
        if (Def_DWM_mixL > Def_DWM * 0.95)==1:   # originally 0.95
            #print '(Def_DWM_mixL > Def_DWM * 0.95)==1'
            break
        elif (b_loop == (meta.dr_mixl * meta.lr_mixl) -1)==1:
            # if this condition is checked, it supposed that the wake is more wider than mixL domain
            # so we should increase lr_mixL
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print 'you should increase lr_mixl in meta. '
            print 'the wake is wider than mixL domain'
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            break
    # if j==2:
    #     elapsed_solve = time.time() - tsolve
    #     print 'Solver matrix computation time %i' % (int(elapsed_solve))
    width[j-1]    =  b_loop
    #print 'width: ', width

    return width


def DWM_eddy_viscosity(mfor,meta,width,visc_wake,visc_wake1,visc_wake2,F1_vector,F2_vector,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,j):
    """Function that calculate the eddy viscosity

        Inputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        mfor (instance of class): Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
        F1_vector (np.array float): DWM filter functions F 1 governing the development of turbulent stresses ref[6-7], eq 3 of [5]
        F2_vector (np.array float): DWM filter functions F 2 governing the development of turbulent stresses ref[6-7], eq 3 of [5]
        visc_wake1 (np.array float): initialized contribution due to ambient turbulence for the eddy viscosity ref[5], Eq 3 first term of right and side
        visc_wake2 (np.array float): initialized contribution due to shear layer of wake deficit for the eddy viscosity ref[5], Eq 3 last term of right and side
        visc_wake (np.array float): eddy viscosity of DWM, combination of visc_wake1 and visc_wake2
        u_star_DEF (np.array float): non dimensional velocity scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
        l_star_DEF (np.array float): non dimensional integral length scale ambient turbulence, which affect the wake deficit evolution (roughly corresponding to eddies smaller than 2D)
        One_div_du_dr_DWW (np.array float): denominator of Eq (5) in ref [5]
        width (np.array int): wake width vector index
        j (int) : index in the main Ainslie forward scheme loop

        Outputs
        ----------
        mfor (instance of class): updated instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
    """


## Calculate eddy viscosity
    #print 'j: ', j
    #print 'meta.lr_mixl*meta.dr_mixl (float): ', meta.lr_mixl*meta.dr_mixl
    meta.nr_mixl = int(meta.lr_mixl*meta.dr_mixl)
    #print 'nr int: ', meta.nr_mixl
    ## Include blend between original Prandtl model and Ainslie to avoid issues when wake turbulence goes to 0.
    ## The largest eddy viscosity at each point is applied.
    # Calculate mean flow gradient - du/dr is created with CDS (apart from 1st and last point)
    mfor.du_dr_DWM[j-1,0]                  = (mfor.U[j-1,1] - mfor.U[j-1,0])/meta.dr
    mfor.du_dr_DWM[j-1,1: meta.nr_mixl-2]  = (mfor.U[j-1,2:(meta.nr_mixl-2)+1] - mfor.U[j-1,0:(meta.nr_mixl-2)-1])/(2*meta.dr)
    mfor.du_dr_DWM[j-1,meta.nr_mixl-1]      = (mfor.U[j-1, meta.nr_mixl-1] - mfor.U[j-1, meta.nr_mixl-2])/meta.dr

    if meta.Keck or meta.previous_sDWM:
        # Blend of mixL and Ainslie eddy visc
        # Keck Visc
        visc_wake1[j - 1, :] = F2_vector[j - 1] * meta.k2 * (meta.vr_mixl[width[j - 1] - 1] / meta.R_WTG) ** 2 * np.abs(
            mfor.du_dr_DWM[j - 1, :])
        # Ainslie Visc
        visc_wake2[j - 1, :] = F2_vector[j - 1] * meta.k2 * (meta.vr_mixl[width[j - 1] - 1] / meta.R_WTG) * (
        1.0 - np.min(mfor.U[j - 1, :]))
        # Take the max of two
        visc_wake[j - 1, :] = np.maximum(visc_wake1[j - 1, :], visc_wake2[j - 1, :])
        # max operator is included in the eddy viscosity formulation to avoid underestimating the
        # turbulent stresses at locations where the velocity gradient of the deficit du_dr approaches zero
        # Atmospheric eddy visc as u*l*, yields total eddy viscosity

        visc_norm_factor      = 6.3918 # Applied to use Keck et al. [3] calibration
        mfor.visc[j-1,:]           = F1_vector[j-1]*meta.k1*visc_norm_factor*u_star_DEF*l_star_DEF + visc_wake[j-1,:]

    if meta.Madsen:
        # Madsen visc presented November 2010 'Calibration and Validation of DWM for implementation in an aero code'
        if not meta.without_filter_functions:
            visc_wake[j - 1, :] = F2_vector[j-1]*meta.k2_Madsen * (meta.vr_mixl[width[j - 1] - 1] / meta.R_WTG) * (1.0 - np.min(mfor.U[j - 1, :]))
            mfor.visc[j-1,:] = F1_vector[j-1] *meta.k_amb_Madsen * meta.mean_TI_DWM + visc_wake[j - 1,:]

        # Madsen visc presented November 2008 'Wake deficit and turbulence simulated with two models...'
        elif meta.without_filter_functions:
            visc_wake[j - 1, :] = meta.k2_Madsen * (meta.vr_mixl[width[j - 1] - 1] / meta.R_WTG) * (
            1.0 - np.min(mfor.U[j - 1, :]))
            mfor.visc[j - 1, :] = meta.k_amb_Madsen * meta.mean_TI_DWM + visc_wake[j - 1, :]
        #raise Exception('Eddy Visc Not Implemented for now')

    if meta.Larsen:
        visc_wake[j - 1, :] = F2_vector[j - 1] * meta.k2_Madsen * (meta.vr_mixl[width[j - 1] - 1] / meta.R_WTG) * (1.0 - np.min(mfor.U[j - 1, :]))
        mfor.visc[j - 1, :] = F1_vector[j - 1] * meta.F_amb * meta.k_amb_Larsen * meta.mean_TI_DWM + visc_wake[j - 1, :]


    if meta.Keck or meta.previous_sDWM:
        ## Include contribution from atmospheric boundary layer on DWM
        ##  turbulent stresses. This effect is taken into account by:
        # 1. Calculate the azimuthally averaged local gradient (du/dr tot) acting of the eddy viscosity as a combination of du/dr in the DWM model and du/dz from ABL
        # 2. The du/dr contribution is constant in azimuthal direction. The du/dz part is assumed linear, which gives a sinus curve in a du/dr system
        # 3. Instead of manipulating the velocity field, the eddy viscosity is increased by a "du/dr_total / du/dr_DWM"
        #=> Visc*        =  Visc * du/dr_total / du/dr_DWM
        #   => Turb_Stress  =  Visc* * du/dr_DWM = Visc * (du/dr_total / du/dr_DWM) * du/dr_DWM =  Visc * du/dr_total
        # 4. "Wiener filter" is used to avoid problems when du/dr = 0, idea:  1/f(x) ~= f(x) / (f(x)^2 + k)

        # Calculate total mean flow gradient - adds shear contribution via
        # sinus function. This gets the stresses right, but sign is wrong in
        #regions where du/dr_DWM - sign of du/dz_ABL is negative
        # notations as per Keck et al. [3].
        du_dr_DWM_du_dz=np.array(np.absolute(mfor.du_dr_DWM[j-1,:]) / mfor.Shear_add_du_dz ,dtype=complex)
        alfa_1      = np.arcsin(du_dr_DWM_du_dz)
        alfa_2      = pi - alfa_1
        alfa_1=np.asarray([abs(x) for x in alfa_1])
        alfa_2=np.asarray([abs(x) for x in alfa_2])
        mfor.du_dr_tot[j-1,0: meta.nr_mixl] = ( np.absolute(mfor.du_dr_DWM[j-1,:]) *2.0*pi +\
        ((np.absolute(mfor.du_dr_DWM[j-1,:]) < mfor.Shear_add_du_dz) * 2.0 * \
        (mfor.Shear_add_du_dz*2.0*np.cos(alfa_1) - np.absolute(mfor.du_dr_DWM[j-1,:])*(alfa_2 - alfa_1) ) ) ) / (2.0*pi)
        # du/dr_DWM block of area
        # condition for added shear gradient (if du/dr_DWM >= du/dz_ABL there are no contribution)
        # Area A1 + A2 in figure XXX
        # Scaling from area to gradient
        k_wiener                 = 2.0*mfor.Shear_add_du_dz * meta.dr**2
        One_div_du_dr_DWM[j-1,:] = mfor.du_dr_DWM[j-1,:] / (mfor.du_dr_DWM[j-1,:]**2 + k_wiener)
        visc_fac                 = np.maximum(1.0, (mfor.du_dr_tot[j-1,:] * np.fabs(One_div_du_dr_DWM[j-1,:])))
        mfor.visc[j-1,:]              = mfor.visc[j-1,:] * visc_fac

    return mfor


def DWM_velocity_solver(mfor,meta,j):
    """Function that defines tridiagonal matrix to solve the NS equations Eq 1 and in Ref [5]
       The momentum equation is discretized using a second order central difference scheme in radial direction and a
       first order upwind scheme in flow direction

        (1) Solve the momentum equation for the streamwise velocity component at all radial positions explicitly,
        by using the value of the radial velocity component and the eddy viscosity from the previous location upstream.
        This yields a tri- diagonal equation system where all the coefficients are known, which can easily
        be solved by any tridiagonal ma- trix algorithm.
        (2) Once the streamwise velocity is known, the radial velocity for all radial positions can be updated using the
        continuity equation
        (3) The eddy viscosity for all radial positions is updated using Eq. (3) in Ref [6]
        (4) March to the next downstream location and repeat steps 1-3.

        Inputs
        ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        mfor (instance of class): Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model
        j (int) : index in the main Ainslie forward scheme loop

        Outputs
        ----------
        HL (np.array float)
        mat (np.array float)

    """

    ####### SAME IN VECTORIZED FORMAT FAST !!!
    # vectorization of streamwise loop
    # t = time.time()
    ind_R=range(1,int(meta.lr_mixl*meta.dr_mixl-1),1)
    ind_R_p=[z+1 for z in ind_R]
    ind_R_m=[z-1 for z in ind_R]

    HL = np.zeros((meta.nr_mixl))
    mat = np.zeros(((meta.nr_mixl),(meta.nr_mixl)))

    # Input BC for wake center
    HL[0]     = (mfor.U[j-1,0]**2     / meta.dz_mixl)
    HL[ind_R] = (mfor.U[j-1,ind_R]**2/meta.dz_mixl)
    HL[meta.nr_mixl-1]       = (mfor.U[j-1,meta.nr_mixl-1]/ meta.dz_mixl)


    mat[0,0]  =  mfor.U[j-1,0]/meta.dz_mixl     + (2.0*mfor.visc[j-1,0]/(meta.dr**2))
    mat[1,0]  = -(2.0*mfor.visc[j-1,0] /(meta.dr**2))

    VL11=np.zeros((meta.nr_mixl)-1)
    VL21=np.zeros((meta.nr_mixl)-1)
    VL31=np.zeros((meta.nr_mixl)-1)
    VL41=np.zeros((meta.nr_mixl)-1)
    VL12=np.zeros((meta.nr_mixl)-1)
    VL13=np.zeros((meta.nr_mixl)-1)
    VL22=np.zeros((meta.nr_mixl)-1)
    VL23=np.zeros((meta.nr_mixl)-1)
    VL33=np.zeros((meta.nr_mixl)-1)
    VL43=np.zeros((meta.nr_mixl)-1)
    VL1=np.zeros((meta.nr_mixl)-1)
    VL2=np.zeros((meta.nr_mixl)-1)
    VL3=np.zeros((meta.nr_mixl)-1)

    # Calculation of U for the wake body
    # mfor.visc[j-1,range(0,int(meta.lr_mixl*meta.dr_mixl),1)+1]
    VL11[ind_R]             = -mfor.V[j-1,ind_R]      / (2.0*meta.dr)
    VL21[ind_R]             = mfor.visc[j-1,ind_R]    / (2.0*meta.vr_mixl[ind_R]*meta.dr)
    VL31[ind_R]             = -mfor.visc[j-1,ind_R]   / (meta.dr**2)
    VL41[ind_R]             = (mfor.visc[j-1,ind_R_p] - mfor.visc[j-1,ind_R_m])  / (2*meta.dr)**2 # new term due to d(nu_t)/dr dependence
    VL12[ind_R]             = mfor.U[j-1,ind_R]       / (meta.dz_mixl)
    VL22[ind_R]             = +2.0*mfor.visc[j-1,ind_R] / (meta.dr**2)
    VL13[ind_R]             = mfor.V[j-1,ind_R]       / (2.0*meta.dr)
    VL23[ind_R]             = -mfor.visc[j-1,ind_R]   / (2.0*meta.vr_mixl[ind_R]*meta.dr)
    VL33[ind_R]             = -mfor.visc[j-1,ind_R]   / (meta.dr**2)
    VL43[ind_R]             = -(mfor.visc[j-1,ind_R_p] - mfor.visc[j-1,ind_R_m])  / (2.0*meta.dr)**2 # new term due to d(nu_t)/dr dependence
    VL1[ind_R]              = VL11[ind_R] + VL21[ind_R] + VL31[ind_R] + VL41[ind_R]
    VL2[ind_R]              = VL12[ind_R] + VL22[ind_R]
    VL3[ind_R]              = VL13[ind_R] + VL23[ind_R] + VL33[ind_R] + VL43[ind_R]

    # build the matrix for X =A/B
    mat[ind_R_m,ind_R] = VL1[ind_R]
    mat[ind_R ,ind_R] = VL2[ind_R]
    mat[ind_R_p,ind_R] = VL3[ind_R]

    # Input BC for wake edge
    VL1                   = 0.0
    VL2                   = 1.0/meta.dz_mixl
    mat[meta.nr_mixl-2, meta.nr_mixl-1]  = VL1
    mat[meta.nr_mixl-1, meta.nr_mixl-1]  = VL2

    mat=mat.T
    HL=HL.T
    # elapsed = time.time() - t
    # print 'New loop', 100.*elapsed
    return HL, mat, HL, mat


def DWM_calc_mixL(meta,aero,mfor):
    """ Main Ainslie - mixing length (Keck et al) function that compute the wake deficit as function of downstream
    distance in the meandering frame of reference

    (1) Solve the momentum equation for the streamwise velocity component at all radial positions explicitly, by using
    the value of the radial velocity component and the eddy viscosity from the previous location upstream.
    This yields a tri- diagonal equation system where all the coefficients are known, which can easily be solved by any
    tridiagonal ma- trix algorithm.
    (2) Once the streamwise velocity is known, the radial velocity for all radial positions can be updated using the
    continuity equation
    (3) The eddy viscosity for all radial positions is updated using Eq. (3) in Ref [6]
    (4) March to the next downstream location and repeat steps 1-3

    Inputs
    ----------
        meta (instance of class): Instance of class Meta holding DWM core variables
        aero (instance of class): Instance of class Aero holding BEM-aero core variables
        mfor (instance of class): Instance of class Mfor holding the meandering frame of reference scalars used by the Ainslie model

    Outputs
    ----------
        mfor (instance of class): Updated instance of class Mfor holding the velocity deficit

    """
    print '# -------------------------- # mixL DEFICIT MODULE PROCESSING # ------------------------------------------ #'
    print '--------------------------  Eddy viscosity model used: ' + meta.EV_model + ' -----------------------------'
    save_wakewidth_data = False
    i_p = 0
    if save_wakewidth_data:
        print "Saving wake width data"
        wake_width = []
#    t = time.time()
    # Init arrays and constants for mfor_mixL module
    meta,mfor,F1_vector,F2_vector,visc_wake1,visc_wake2,visc_wake,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,width=DWM_init_calc_mixl(meta,aero,mfor)

    #  Start x-stepping & solving of BLE equations
    # dr_DWM        = 1.0/meta.dr_mixl
    # r_vec_DWM     = np.arange(dr_DWM/2.0,meta.lr_mixl-dr_DWM/2.0+dr_DWM,dr_DWM)
    # dA_DWM        = np.array([pi*r_vec_DWM[1:]**2-pi*r_vec_DWM[0:-1]**2])
    b_loop=meta.dr_mixl
    for j in np.arange(1,len(meta.vz_mixl),1):
        ## Calculating wake width
        width=DWM_calc_wake_width(mfor,width,b_loop,meta,j)

        ## Calculate eddy viscosity
        mfor=DWM_eddy_viscosity(mfor,meta,width,visc_wake,visc_wake1,visc_wake2,F1_vector,F2_vector,u_star_DEF,l_star_DEF,One_div_du_dr_DWM,j)

        ## Calculation matrix
        HL, mat, HL, mat=DWM_velocity_solver(mfor,meta,j)

        ## Solve for U
        mfor.U[j,:] = linalg.solve(mat, HL)

        for i in np.arange(0,meta.nr_mixl-1,1):
            mfor.V[j,i+1] = (meta.vr_mixl[i] / meta.vr_mixl[i+1]) * mfor.V[j,i] - (meta.dr/(2*meta.dz_mixl))*( (mfor.U[j,i+1] - mfor.U[j-1,i+1]) + \
            (meta.vr_mixl[i] / meta.vr_mixl[i+1])*((mfor.U[j,i] - mfor.U[j-1,i])) )

        if meta.Keck:
            #--------------------------------------------------------------------------------------------------------------#
            # CHECK THIS PART FOR MADSEN, LARSEN
            # POST PROCESSING SIGNAL: Turbulent stress
            mfor.Turb_Stress_DWM[j-1,:]           = mfor.visc[j-1,:]  * mfor.du_dr_DWM[j-1,:]


            # POST PROCESSING SIGNAL: TI_DWM formulated based on the relation derived between u'v' and u'u' (see Keck et al. [3])
            x_uw_wake      = 1.0
            C_uw_wake      = (0.7550 - meta.mean_TI_DWM*1.75) / 2 # Article states: "C_uw_wake = 0.3", but this seems to give better results (=> 0.3 for TI=8.6%)
            mfor.TI_DWM[j-1,:]  = np.sqrt( np.abs( (1.0 / (x_uw_wake * C_uw_wake)) * mfor.Turb_Stress_DWM[j-1,:]) )


        mfor.WakeW = meta.vr_mixl[width-1]

        # Plot deficit and derivative for Madsen Approach
        if meta.WaT_detail and meta.vz_mixl[j-1]< 3.3*2 < meta.vz_mixl[j]:

            # MADSEN Scaling for Turbulence:
            km1 = 0.6
            km2 = 0.035
            Udef = mfor.U[j-1,:]
            # derive U by r:
            dUdef_dr_r = [0] + [(Udef[i_r + 1] - Udef[i_r-1]) / (meta.vr_mixl[i_r + 1] - meta.vr_mixl[i_r-1]) for i_r
                             in range(1,len(meta.vr_mixl) - 1)] + [0]
            #dUdef_dr_r = [0] + [(1.5*Udef[i_r + 1] - 2.*Udef[i_r]+0.5*Udef[i_r-1]) / (meta.vr_mixl[i_r + 1] - meta.vr_mixl[i_r - 1]) for
             #                   i_r
              #                  in range(1, len(meta.vr_mixl) - 1)] + [0]
            #dU_dr_r = np.array(dUdef_dr_r)
            deficit_depth = 1-Udef
            derivative_deficit = dUdef_dr_r
            km_r = np.abs(1-Udef) * km1 + np.abs(dUdef_dr_r) * km2

            plt.figure('derivative and deficit at 3.3D')
            plt.title('at '+str(meta.vz_mixl[j-1])+ ' [R]')
            plt.plot(deficit_depth, meta.vr_mixl, label='Deficit Depth')
            plt.plot(derivative_deficit, meta.vr_mixl, label='dUdef')
            plt.plot(Udef, meta.vr_mixl, label='Udef')
            plt.legend()
            plt.show()
    #print 'WakeWidth shape', wake_width.shape
        #mfor.width = width
    # print mfor.debugV[:,:]-mfor.V[22,:]

    if save_wakewidth_data:
        print np.hstack((np.array(meta.vz_mixl[1:]).reshape(len(meta.vz_mixl[:-1]),1), mfor.WakeW))
        np.save('C:/Users/augus/Documents/Stage/Codes/Mann_Turbulence/Result/Wake_Expansion/Ainslie_Wake_Expansion', np.hstack((np.array(meta.vz_mixl[:-1]).reshape(len(meta.vz_mixl[:-1]),1), mfor.WakeW)))
        print 'wake width data saved... '
#    elapsed = time.time() - t
#    print 'Velocity model computation time %i' % (int(elapsed))
    if meta.Deficit_Process_Detail:
        reverse_axis = True

        # ------------------------------ # WS plot # ------------------------------ #
        plt.figure('Axial Velocity Output from mixL domain (MFOR)')
        plt.title('Axial Velocity Output from mixL domain (MFOR)')
        #plt.xlim(1.,0.)
        plt.ylabel('vr (polar discretization)[R]'), plt.xlabel('U [U0]')
        for i_z in np.arange(0, meta.nz, 1):
            plt.plot(mfor.U[meta.vz[i_z], :], meta.vr_mixl, label='at Turbine '+ str(7-i_z))

        plt.legend(), plt.show()

        plt.figure(2)
        plt.title('WS at different downstream distances')
        L_observer = [1,3,6,10, 16,32,48,64] # in R, strictly ascending location
        i_p = 0
        for j in np.arange(1, len(meta.vz_mixl), 1):
            if meta.vz_mixl[j - 1] <= L_observer[i_p] < meta.vz_mixl[j]:
                i_p = i_p + 1
                Udef = mfor.U[j-1,:]

                if reverse_axis:
                    #plt.plot(Udef, meta.vr_mixl, label = str(meta.vz_mixl[j-1])[:3]+'-'+str(meta.vz_mixl[j])[:3]+' [R]')
                    plt.plot(Udef, meta.vr_mixl, label='at '+str(meta.vz_mixl[j - 1])[:3]+' [R]')
                else:
                    plt.plot(meta.vr_mixl, Udef,label=str(meta.vz_mixl[j - 1])[:3] + '-' + str(meta.vz_mixl[j])[:3] + ' [R]')
                if i_p == len(L_observer):
                    break
        if reverse_axis:
            plt.xlabel('[U0]'), plt.ylabel('[R]')
        else:
            plt.ylabel('[U0]'), plt.xlabel('[R]')
        plt.legend()
        plt.show()

        # ------------------------------- # TI plot # -------------------------- #
        # TI at each Turbine Location
        if meta.Keck:
            plt.figure('Axial TI Output from mixL domain (MFOR)')
            plt.title('Axial TI Output from mixL domain (MFOR)')
            plt.ylabel('vr (polar discretization)[R]')
            plt.xlabel('TI')
            for i_z in np.arange(0, meta.nz, 1):
                plt.plot(mfor.TI_DWM[meta.vz[i_z], :], meta.vr_mixl, label='at Turbine '+ str(7-i_z))

            plt.legend(), plt.show()

        # TI at some downstream location
            plt.figure(3)
            plt.title('TI at different downstream distances')
            L_observer = [1, 3, 10, 16, 32, 48, 64]  # in R, strictly ascending location
            i_p = 0
            for j in np.arange(1, len(meta.vz_mixl), 1):
                if meta.vz_mixl[j - 1] <= L_observer[i_p] < meta.vz_mixl[j]:
                    i_p = i_p + 1
                    TI = mfor.TI_DWM[j - 1, :]
                    if reverse_axis:
                        #plt.plot(TI, meta.vr_mixl, label = str(meta.vz_mixl[j-1])[:3]+'-'+str(meta.vz_mixl[j])[:3]+' [R]')
                        plt.plot(TI, meta.vr_mixl,
                                 label='at '+str(meta.vz_mixl[j - 1])[:3] + '[R]')

                    else:
                        plt.plot(meta.vr_mixl, TI, label=str(meta.vz_mixl[j - 1])[:3] + '-' + str(meta.vz_mixl[j])[:3] + ' [R]')
                    if i_p == len(L_observer):
                        break
            TI_amb = [meta.mean_TI_DWM for i in meta.vr_mixl]
            if reverse_axis:
                plt.plot(TI_amb, meta.vr_mixl, label='TI_Amb')
                plt.xlabel('TI []'), plt.ylabel('[R]')
            else:
                plt.plot(meta.vr_mixl, TI_amb, label = 'TI_Amb')
                plt.ylabel('TI []'), plt.xlabel('[R]')
            plt.legend()
            plt.show()

        if meta.AINSLIE_EV_details:
            reverse_axis = False
            print 'Plot: We restrict the x abscisse to 2.5'
            plt.figure()
            plt.title('Eddy Viscosity at different downstream distances')
            L_observer = [1, 3, 10, 16, 32, 64]  # in R, strictly ascending location
            i_p = 0
            for j in np.arange(1, len(meta.vz_mixl), 1):
                if meta.vz_mixl[j - 1] <= L_observer[i_p] < meta.vz_mixl[j]:
                    i_p = i_p + 1
                    visc = mfor.visc[j, :]
                    if reverse_axis:
                        #plt.plot(visc, meta.vr_mixl,label=str(meta.vz_mixl[j - 1])[:3] + '-' + str(meta.vz_mixl[j])[:3] + ' [R]')
                        plt.plot(visc, meta.vr_mixl, label='at '+str(meta.vz_mixl[j - 1])[:3] + ' [R]')
                    else:
                        plt.plot(meta.vr_mixl, visc,label=str(meta.vz_mixl[j - 1])[:3] + '-' + str(meta.vz_mixl[j])[:3] + ' [R]')
                    if i_p == len(L_observer):
                        break

            if reverse_axis:
                plt.ylabel('Radial position [R]'), plt.xlabel('eddy viscosity [-]')
            else:
                plt.xlabel('Radial position [R]'), plt.ylabel('eddy viscosity [-]')
    if meta.WaT_detail:
        reverse_axis = False
        print 'Plot: We restrict the x abscisse to 2.5'
        plt.figure()
        plt.title('Radial scaling Factor at different downstream distances')
        L_observer = [1, 3, 10, 16, 32, 64]  # in R, strictly ascending location
        i_p = 0
        for j in np.arange(1, len(meta.vz_mixl), 1):
            if meta.vz_mixl[j - 1] <= L_observer[i_p] < meta.vz_mixl[j]:
                i_p = i_p + 1
                if meta.Madsen or meta.Larsen:
                    Udef = mfor.U[j-1,:]
                    # derive U by r:
                    dUdef_dr_r = [0] + [(Udef[i_r + 1] - Udef[i_r - 1]) / (meta.vr_mixl[i_r + 1] - meta.vr_mixl[i_r - 1])
                                        for i_r in range(1, len(meta.vr_mixl) - 1)] + [0]

                    deficit_depth = 1 - Udef
                    derivative_deficit = dUdef_dr_r

                    kmt_r = np.abs(deficit_depth) * meta.km1 + np.abs(derivative_deficit) * meta.km2

                if meta.Keck:
                    kmt_r = mfor.TI_DWM[j - 1, :]/meta.mean_TI_DWM
                plt.plot(kmt_r, meta.vr_mixl,
                         label=str(meta.vz_mixl[j - 1])[:3] + '-' + str(meta.vz_mixl[j])[:3] + ' [R]')
                if i_p == len(L_observer):
                    break
        plt.ylabel('Radial position [R]'), plt.xlabel('[-]')
        #plt.xlim(0., 2.5)
        plt.legend(), plt.show()
    print '# -------------------------- # mixL DEFICIT MODULE PROCESS ENDED # --------------------------------------- #'
    return mfor