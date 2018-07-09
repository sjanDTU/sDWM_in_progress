# -*- coding: utf-8 -*-
"""
Created on Tue May  5 08:44:05 2015

@author: Ewan Machefaux
# Based on the BEM of Dr. Emmanuel Branlard
"""

# import matplotlib as plt
import time
from sys import exit
import numpy as np
from math import cos, atan2, acos, exp, sqrt, sin
from cBEM import *
from os import chdir, listdir, path, getcwd
from glob import glob
from DWM_misc import find_intersections, unique_rows, unique2d
import os
# import matlab.engine
import matplotlib.pyplot as plt
import matplotlib._cntr as cntr
from scipy.optimize import fmin_cobyla

def getOptimalLambdaPitch(meta,U0,derating,Spec,Rotor):

    Pd=Spec.P_rated *derating
    # Rotor.SweptArea = pi * (Rotor.R * cos(Rotor.cone * pi / 180.)) ** 2
    # Rotor.SweptArea = pi * (Rotor.R * cos(Rotor.cone * pi / 180.)) ** 2
    # Rotor.R_coned = Rotor.R * cos(Rotor.cone * pi / 180.)  # Rotor coned radius [m]
    CPde=(Pd*1000)/(0.5*1.225*U0**3*(Rotor.SweptArea))
    maxCP=np.amax(np.amax(meta.WTG_spec.CP3))
    Pref=meta.WTG_spec.get_P(U0)*1000
    CPrefe=Pref/(0.5*1.225*((U0)**3)*(Rotor.SweptArea))
    CPref=np.minimum(CPrefe,0.999995*maxCP)
    CPd=np.minimum(CPref,CPde)
    CPde=CPd
    # print 'Cpd is', CPde
    # print 'Max cp is', maxCP
    # print 'Pd is', Pd


    #########################################
    # Normal derating following the existing pitch/lambda curve
    #########################################
    flag=False
    lambdad_ref=(pi*Spec.RPM*Rotor.R)/(30*Spec.WS)
    segs_lambda=np.zeros((len(Spec.Pitch),2),dtype=float)
    segs_lambda[:,0]=Spec.Pitch
    segs_lambda[:,1]=lambdad_ref
    segs_lambda_unique=segs_lambda

    # method 1
    # segs_lambda_unique=unique2d(segs_lambda_unique)
    # segs_lambda_unique=unique_rows(segs_lambda_unique)
    # segs_lambda_unique=segs_lambda_unique[segs_lambda_unique[:,1].argsort()[::-1]]


    # method 2
    b = np.ascontiguousarray(segs_lambda_unique).view(np.dtype((np.void, segs_lambda_unique.dtype.itemsize * segs_lambda_unique.shape[1])))
    _, idx, inv = np.unique(b, return_index=True, return_inverse=True)
    unique_a = segs_lambda_unique[idx]
    segs_lambda_unique=unique_a[unique_a[:, 1].argsort()[::-1]]
    # segs_lambda_unique=segs_lambda_unique_if


    # This is an add-on for the NY2 turbine, because the intersection algo does not work with straight vertical line in the below rated region at constant pitch
    segs_lambda_unique_lin=segs_lambda_unique[[2,13],:]
    segs_lambda_unique_lin[0]=segs_lambda_unique_lin[0]*1.00000001
    segs_lambda_unique_lin[1]=segs_lambda_unique_lin[1]*0.99999999
    # print 'segs_lambda_unique_lin', segs_lambda_unique_lin
    # print 'segs_lambda_unique', segs_lambda_unique

    Pref=meta.WTG_spec.get_P(U0)*1000
    CPrefe=Pref/(0.5*1.225*((U0)**3)*(Rotor.SweptArea))
    CPref=np.minimum(CPrefe,maxCP)
    CPd=np.minimum(CPref,CPde)
    # print 'CPref %4.10f / Cpd %4.10f / WS %4.2f / de %4.2f'% (CPref,CPde,U0,derating)
    l_d = CPd
    nlist = meta.WTG_spec.CPc.trace(l_d, l_d, 0)
    segs_d = nlist[:len(nlist)//2]

    # segs_lambda_unique=unique_rows(segs_lambda)
    # segs_lambda_unique=unique_rows(segs_lambda)
    # segs_lambda_unique=segs_lambda_unique[segs_lambda_unique[:,1].argsort()[::-1]]
    try:
        # segs_d = nlist[:len(nlist)//2]
        len(nlist[0])
    except:
        print 'Exception was raised, no interesection between CPd and lambda_ref, probably HIGH CP requested'
        theta_p=np.interp(U0,Spec.WS,Spec.Pitch)
        TSR=np.interp(U0,Spec.WS,Spec.RPM*pi*Rotor.R)/(30*U0)
        # TSR=np.amax(lambdad_ref)
        # print lambdad_ref
        # print TSR
        # print np.argmax(lambdad_ref)
        # theta_p=Spec.RPM[np.argmax(lambdad_ref),2]
        # TSR=lambdad_ref[0]
        # theta_p=Spec.Pitch[0]
        flag=True
    if flag is not True:
        try:
            x, y = find_intersections(segs_d[0],segs_lambda_unique)
            if len(x)==0:
                # print 'Couldn''t find interesection, use the linear approach'
                raise ValueError('Couldn''t find interesection, use the linear approach')
            else:
                # print x,y
                # plt.figure()
                # plt.plot(segs_lambda_unique[:,0],segs_lambda_unique[:,1])
                # plt.plot(segs_d[0][:,0],segs_d[0][:,1],'r-')
                # plt.plot(x,y,'go')
                theta_p=x[0]
                TSR= y[0]
        except:
            # print 'Exception was raised, no interesection between CPd and lambda_ref, probably low CP'
            x, y = find_intersections(segs_d[0],segs_lambda_unique_lin)
            # print x,y
            # plt.figure()
            # plt.plot(segs_lambda_unique_lin[:,0],segs_lambda_unique_lin[:,1])
            # plt.plot(segs_d[0][:,0],segs_d[0][:,1],'r-')
            # plt.plot(x,y,'go')
            theta_p=x[0]
            TSR= y[0]

    try:
        x, y = find_intersections(segs_lambda,segs_d[0])
        theta_p=np.amax(x)
        TSR= y[np.argmax(x)]
    except:
        print 'Exception was raised, no interesection between CPd and lambda_ref'
        # theta_p=np.nan
        # TSR=np.nan
        theta_p=np.interp(U0,Spec.WS,Spec.Pitch)
        TSR=np.interp(U0,Spec.WS,Spec.RPM*pi*Rotor.R)/(30*U0)
    #########################################
    #########################################
    # RPM=(TSR*30.*U0)/(pi*Rotor.R)
    # Pitch=theta_p


    #########################################
    # Old method where the tangent is not accurately found and where an actual contour plot of the CT was done (memory leak on ipython notebook)
    #########################################
    # lambdad=[]
    # pitchd=[]
    # for iCt in range(0,len(meta.WTG_spec.CTc.levels)):
    #     datCT= meta.WTG_spec.CTc.collections[iCt].get_paths()[0].vertices
    #     x, y = find_intersections(datCT,segs_d[0])
    #     if len(x) == 0:
    #         pass
    #     else:
    # #         print x, y
    #         lambdad.append(x)
    #         pitchd.append(y)
    #         break
    # strategy='low omega' # low pitch
    # # strategy='low pitch'
    # if strategy == 'low omega':
    #     try:
    #         TSR=min(lambdad[0])
    #         theta_p= pitchd[0][np.argmin(lambdad[0])]
    #     except:
    #         print 'cannot find optimum set points, check derating factor maybe its too low.. in the mean time is use the nominal pitch and rpm'
    #         print 'Cpd is', CPd
    #         RPM=np.interp(U0,Spec.WS,Spec.RPM)
    #         TSR= (RPM*pi*Rotor.R)/(30.*U0)
    #         theta_p=np.interp(U0,Spec.WS,Spec.Pitch)
    # elif strategy == 'low pitch':
    #     try:
    #         theta_p=min(pitchd[0])
    #         TSR= lambdad[0][np.argmin(pitchd[0])]
    #     except:
    #         print 'cannot find optimum set points, check derating factor maybe its too low.. in the mean time is use the nominal pitch and rpm'
    #         print 'Cpd is', CPd
    #         RPM=np.interp(U0,Spec.WS,Spec.RPM)
    #         TSR= (RPM*pi*Rotor.R)/(30.*U0)
    #         theta_p=np.interp(U0,Spec.WS,Spec.Pitch)
    #########################################
    #########################################

    #########################################
    # Better tangente to the CT and CP curve
    #########################################
    # lambdad=[]
    # pitchd=[]
    # nlist = meta.WTG_spec.CPc.trace(CPde, CPde, 0)
    # segs_d = nlist[:len(nlist)//2]
    # # print segs_d
    # for iCt in np.arange(0.7*CPde,1.,0.0001): # change lower and upper bounds for speed based on CT, CP curve (implemented soon)
    #     nlist=meta.WTG_spec.CTc.trace(iCt,iCt,0)
    #     segs_ct = nlist[:len(nlist)//2]
    # #     datCT= CTc.collections[iCt].get_paths()[0].vertices
    #     x, y = find_intersections(segs_ct[0],segs_d[0])
    #     if len(x) == 0:
    #         pass
    #     else:
    # #         print x, y
    #         pitchd.append(x)
    #         lambdad.append(y)
    #         # print iCt
    #         break
    # # print lambdad
    # # print pitchd
    # TSR=np.mean(lambdad)
    # theta_p=np.mean(pitchd)
    #########################################
    #########################################

    RPM=(TSR*30.*U0)/(pi*Rotor.R)
    Pitch=theta_p
    # print 'New RPM', RPM
    # print 'New pitch', Pitch
    # print 'Standard RPM', np.interp(U0,Spec.WS,Spec.RPM)
    # print 'Standard Pitch', np.interp(U0,Spec.WS,Spec.Pitch)
    return RPM, Pitch

#########################################################################################
#########################################################################################
def getInduction(ngrid, sWT, Format, U0, meta, ID_waked, **kwargs):
    opt = 'python'
    # opt = 'matlab'
    derating = kwargs.get('derating')
    # Initialization
    PcDataAll, AeData, Rotor, Env, WT, Spec, State, Algo, Misc = fInitWT(sWT, Format, '../WT-data/')
    Rotor = fSetRotorGrid(ngrid, Rotor, Env, WT, Spec, State, Algo, Misc)

    # if derating == 1.:
    RPM=np.interp(U0,Spec.WS,Spec.RPM)
    Pitch=np.interp(U0,Spec.WS,Spec.Pitch)
    # RPM, Pitch= getOptimalLambdaPitch(meta,U0,derating,Spec,Rotor)
    # float("{0:.2f}".format(meta.WTG_spec.get_P(meta.mean_WS_DWM)))
    Sim = InitSim()
    parSim = {'Name': sWT,
          'WT': sWT,
          'WS': U0,
          'RPM': RPM,
          'PITCH': Pitch,
          'RPM_opt' : np.interp(U0,Spec.WS,Spec.RPM),
          'PITCH_opt' : np.interp(U0,Spec.WS,Spec.Pitch)}
    Sim.Set(**parSim)
    # print Spec.RPM

    Wind = InitWind()
    parWind = {'V0': np.array([0., 0., U0]),
           'WT': sWT,
           'WS': U0,
           'RPM': RPM,
           'PITCH': Pitch}
    Wind.Set(**parWind)
    # print Spec.P_rated
    # print 'U0 is', U0
    # print 'RPM is', RPM
    # print 'Pitch is', Pitch
    Algo = InitAlgo()
    BEM = fBEMsteady(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc)
    # print 'float("{0:.2f}".format(meta.WTG_spec.P_rated *derating))',float("{0:.2f}".format(meta.WTG_spec.P_rated *derating))
    # print 'float("{0:.2f}".format(BEM.Power/1000.)', float("{0:.2f}".format(BEM.Power/1000.))
    if float("{0:.2f}".format(meta.WTG_spec.P_rated *derating))  >=  float("{0:.2f}".format(BEM.Power/1000.)):
        BEM.derated=False
        pass
        # print 'You have demanded a derating but I produce less anyway'
    elif opt == 'matlab':
        initcd=os.getcwd()
        os.chdir('../BEM/BEM_v05/BEM/v05/Simulations/')
        eng = matlab.engine.start_matlab()
        eng.addpath(r'/home/ewan/Work/Msc/MatlabPath/',nargout=0)

        ## INPUTS TO THE get Optimal routine in matlab
        # derating=0.2774
        # U=9.
        Powcurve=matlab.double(meta.WTG_spec.ref_P.tolist())
        # print Powcurve
        Wscurve= matlab.double(meta.WTG_spec.ref_u.tolist())
        RPMcurve=matlab.double(Spec.RPM.tolist())


        Power,Pitch,RPM,WS,Ud,BEM_a,BEM_r,BEM_CP,BEM_CT =eng.fControllerRef(derating,U0,Wscurve,Powcurve,RPMcurve,nargout=9)
        os.chdir(initcd)
        Sim = InitSim()
        parSim = {'Name': sWT,
              'WT': sWT,
              'WS': U0,
              'RPM': RPM,
              'PITCH': Pitch,
              'RPM_opt' : np.interp(U0,Spec.WS,Spec.RPM),
              'PITCH_opt' : np.interp(U0,Spec.WS,Spec.Pitch)}
        Sim.Set(**parSim)
        ne = Rotor.ne
        BEM = InitBEM(ne)
        BEM.derated=True
        BEM.CP=BEM_CP
        # BEM.CPloc=BEM_CPloc
        BEM.a=np.asarray(BEM_a[0])
        BEM.r=np.asarray(BEM_r[0])
        BEM.R=Rotor.R
        BEM.CT=BEM_CT
        # BEM.RPM = Sim.RPM
        # BEM.PITCH = Sim.PITCH
        BEM.RPM_opt = Sim.RPM_opt
        BEM.PITCH_opt = Sim.PITCH_opt
        BEM.Power=Power
        BEM.PITCH=Pitch
        BEM.RPM=RPM

        print 'Matlab pitch', BEM.PITCH
        print 'Matlab RPM', BEM.RPM
    elif opt == 'python':
        Pd=meta.WTG_spec.P_rated *derating*1000.
        # derated RPM calculation
        iP=np.argmin(abs(meta.WTG_spec.ref_P-meta.WTG_spec.P_rated))
        # Powcurve=meta.WTG_spec.ref_P[0:iP]
        Powcurve=meta.WTG_spec.ref_P
        #print 'iP is', iP
        #print 'Powcurve is', Powcurve
        # WScurve=meta.WTG_spec.ref_u[0:iP]
        WScurve=meta.WTG_spec.ref_u
        # RPMcurve=Spec.RPM[0:iP]
        RPMcurve=Spec.RPM
        Ud   = np.interp(Pd*0.001,Powcurve,WScurve)
        #print 'Ud is', Ud
        RPMd = np.interp(Ud,WScurve,RPMcurve)
        Sim.RPM= RPMd
        #print 'Python Pitch RPMd is', RPMd
        #print 'Pd is ', Pd
        #print 'BEM.Power is', BEM.Power
        x0 =  np.interp(U0/derating,Spec.WS,Spec.Pitch) # init the init derating vector
        # print x0
        # x=fmin_cobyla(objective , x0,cons=[constr1,], args=(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc,Pd),consargs=(derating,), rhobeg=0.05, rhoend=0.000001, iprint=1, maxfun=1000, disp=1, catol=0.002)
        x=fmin_cobyla(objective , x0,cons=[], args=(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc,Pd),consargs=None, rhobeg=0.5, rhoend=0.001, iprint=0, maxfun=1000, disp=0, catol=0.002)
        # x=fFindPitch(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc,Pd)
        # print 'Python Pitch from cobyla is', x
        Sim.PITCH=x
        BEM = fBEMsteady(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc)
        BEM.Ud=Ud
        BEM.derated=True
    ################################################################
    ################################################################
    return BEM
#########################################################################################
#########################################################################################


def objective(x,WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc,Pd):
    Sim.PITCH=x
    BEM=fBEMsteady(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc)
    # print 'BEM.Power-Pd',BEM.Power-Pd
    # print 'SIM.PITCH', Sim.PITCH
    return abs(BEM.Power-Pd)
    # return abs(Ptot-Prated) # maximize power by derating

def constr1(x,derating):
    return x
#
# def constr2(x,*args):
#     return np.min(x - 0.05*np.ones_like(x))

def fBEMsteady(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc):
    cone = Rotor.cone
    # print 'cone is', Rotor.cone
    V0 = Wind.V0[2] * cos(cone * pi / 180.)
    # V0 = Wind.V0[2]
    VHubHeight = Wind.V0[2]
    nB = Rotor.nB
    ne = Rotor.ne
    r = Rotor.r
    # print r
    dr = Rotor.dr
    rhub = Rotor.rhub
    R = Rotor.R
    chord = Rotor.chord
    twist = Rotor.twist
    # print 'raidus is',R
    sigma = chord * nB / (2.0 * pi * r * cos(cone * pi / 180.))
    # Environment
    rho = Sim.rho
    KinVisc = Sim.KinVisc
    pitch = Sim.PITCH
    omega = Sim.RPM * 2. * pi / 60.
    # print omega
    lambda_r = omega * r * cos(cone * pi / 180.) / V0
    lambda_ = omega * R * cos(cone * pi / 180.) / V0
    # algorithm internal paramter
    nbIt = Algo.nbIt
    aTol = Algo.aTol
    Pn = np.zeros((ne))
    Pt = np.zeros((ne))

    BEM = InitBEM(ne)
    ## Main BEM loop
    for e in np.arange(0, ne, 1):
        a = 0.3 * 0.
        aprime = 0.01 * 0.
        Ftip_previous = 1
        for i in np.arange(0, nbIt, 1):
            # --------------------------------------------------------------------------------
            # --- Step 0: Relative wind
            # --------------------------------------------------------------------------------
            # --------------------------------------------------------------------------------
            # ---  Step 1: Wind Components
            # --------------------------------------------------------------------------------
            # if e==0:
            #     print aprime
            Ut = omega * r[e] * (1. + aprime)
            Un = V0 * (1. - a)
            Vrel_norm = np.sqrt(Un ** 2 + Ut ** 2)
            Re = Vrel_norm * chord[e] / KinVisc / 10 ** 6  # Reynolds number in Millions
            # print Ut, Un, Vrel_norm, Re
            # --------------------------------------------------------------------------------
            # --- Step 2: Flow Angle
            # --------------------------------------------------------------------------------
            phi = atan2(Un, Ut) * 180. / pi  # (this returns a positive flow angle) [deg]
            if isinstance(phi, complex):
                print ('Algorithm did not converge : lambda=%.2f beta=%.2f V0=%.2f r=%.2f - it %d' % (
                    lambda_, twist[e], V0, r[e], i))
                phi = 0.
                exit()
            # --------------------------------------------------------------------------------
            # --- Tip loss
            # --------------------------------------------------------------------------------
            Ftip = a * 0. + 1.
            Fperf = a * 0. + 1.
            if Algo.bTipLoss:
                if Algo.TipLossMethod == 'Prandtl':  # originally written Glauert but it's Prandtl
                    if (sin(phi * pi / 180.) > 0.01):
                        # prandtl tip correction
                        Ftip = 2. / pi * acos(exp(-nB / 2. * (R - r[e]) / (r[e] * sin(phi * pi / 180.))))
                elif Algo.TipLossMethod == 'Prandtl....':
                    Ftip = 2. / pi * acos(exp(-nB / 2. * (1. - lambda_r[e] / lambda_) * sqrt(1. + lambda_ ** 2)))
            # here will be implemented Hub losses in the future
            F = Ftip
            # --------------------------------------------------------------------------------
            # --- Step 3: Angle of attack
            # --------------------------------------------------------------------------------
            alpha = phi - (twist[e] + pitch)
            # --------------------------------------------------------------------------------
            # --- Step 4: Profile Data
            # --------------------------------------------------------------------------------
            # print  alpha, twist[e]
            Cl, Cd, Cn, Ct, CnForAI, CtForTI, fs = fAeroCoeffWrap(Rotor, PcDataAll, e, alpha, phi, chord, Vrel_norm, Re,
                                                                  Fperf, WT, Algo, Misc)
            # --------------------------------------------------------------------------------
            # --- Step 5: Induction Coefficients
            # --------------------------------------------------------------------------------
            # a=(V0-Un)/V0;   %induction factor %abs for compatibility with unsteady code
            # Storing last values
            a_last = a
            aprime_last = aprime
            a, aprime, CT_loc = fInductionCoefficients(a_last, [0., 0., Vrel_norm], Un, Ut, [0., 0., V0], [0., 0., V0],
                                                   [0., 0., -a_last * V0], omega, chord[e], F, Ftip, CnForAI, CtForTI,
                                                   lambda_r[e], sigma[e], phi, Algo)

            if (i > 3 and (abs(a - a_last) + abs(aprime - aprime_last)) < aTol):  # used to be on alpha
                break
            Ftip_previous = Ftip
        BEM.nIt[e] = i + 1
        if i == nbIt:
            print('Maximum iterations reached : lambda=%.2f beta=%.2f V0=%.2f r=%.2f' % (lambda_, twist[e], V0, r[e]))
        # --------------------------------------------------------------------------------
        # --- Step 6: Aerodynamic Forces PER LENGTH
        # --------------------------------------------------------------------------------
        # L = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cl
        # D = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cd
        Pn[e] = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Cn
        Pt[e] = 0.5 * rho * Vrel_norm ** 2 * chord[e] * Ct

        BEM.Vrel[e] = Vrel_norm
        BEM.Un[e] = Un
        BEM.Ut[e] = Ut
        BEM.Re[e] = Re
        BEM.F[e] = F
        BEM.Fperf[e] = Fperf
        BEM.a[e] = a
        BEM.a_last[e] = a_last
        BEM.aprime[e] = aprime
        BEM.aprime_last[e] = aprime_last
        BEM.phi[e] = phi
        BEM.alpha[e] = alpha
        BEM.Cl[e] = Cl  # USED TO BE MULTIPLIED BY Fshen, why????
        BEM.Cd[e] = Cd  # SAME
        BEM.Cn[e] = Cn
        BEM.Ct[e] = Ct
        # BEM.CT[e] = CT

    BEM.Pn = Pn
    BEM.Pt = Pt
    BEM.ThrLoc = dr * Pn * cos(cone * pi / 180.)
    BEM.ThrLocLn = Pn * cos(cone * pi / 180.)
    BEM.CTloc = nB * BEM.ThrLoc / (0.5 * rho * VHubHeight ** 2 * (2. * pi * r * cos(cone * pi / 180.) * dr))
    BEM.TqLoc = dr * r * Pt * cos(cone * pi / 180.)
    BEM.TqLocLn = r * Pt * cos(cone * pi / 180.)
    BEM.CQloc = nB * BEM.TqLoc / (
        0.5 * rho * VHubHeight ** 2 * (2. * pi * r * cos(cone * pi / 180.) * dr * r * cos(cone * pi / 180.)))
    BEM.CPloc=BEM.CQloc*lambda_r
    ####  Returning Aerodynamic Forces
    BEM.Torque = nB * getTorqueFromBlade(r, Pt * cos(cone * pi / 180.), R)  # Rotor shaft torque at t in Newton
    BEM.Thrust = nB * getThrustFromBlade(r, Pn * cos(cone * pi / 180.), R)  # Rotor shaft thrust at t in Newton
    BEM.Flap = sum(dr * (Pn * cos(cone * pi / 180.)) * (r - rhub))  # approximate
    BEM.Edge = sum(dr * Pt * (r * cos(cone * pi / 180.)) * (r - rhub))  # approximate
    # print BEM.Un
    # print BEM.Ut
    # print BEM.Torque
    BEM.Power = omega * BEM.Torque
    BEM.CP = BEM.Power / (0.5 * rho * V0 ** 3. * Rotor.SweptArea)
    # print BEM.Thrust
    BEM.CT = BEM.Thrust / (0.5 * rho * V0 ** 2. * Rotor.SweptArea)
    BEM.CQ = BEM.Torque / (0.5 * rho * V0 ** 2. * pi * R ** 3.)
    BEM.Gamma = 0.5 * BEM.Re * BEM.Cl * KinVisc * 10 ** 6.
    BEM.r = r
    BEM.R=R
    BEM.omega=omega
    BEM.uia = V0 * BEM.a
    BEM.uit = omega * r * BEM.aprime
    BEM.RPM = Sim.RPM
    BEM.PITCH = Sim.PITCH
    BEM.RPM_opt = Sim.RPM_opt
    BEM.PITCH_opt = Sim.PITCH_opt
    # print 'BEM.Power is', BEM.Power
    # print 'BEM.CP is', BEM.CP
    # print 'Rotor.SweptArea', Rotor.SweptArea
    # plot the induction factor
    #plt.figure()
    #plt.plot(BEM.r, BEM.a)
    # print BEM.CP
    # print BEM.Power
    # print BEM.RPM
    # print BEM.PITCH
    # print BEM.Vrel
    return BEM
#########################################################################################
#########################################################################################
def getTorqueFromBlade(r0, Pt0, R):
    n = len(r0)
    r = np.zeros((n + 1))
    Pt = np.zeros((n + 1))
    r[0:n] = r0
    r[-1] = R
    Pt[0:n] = Pt0
    Q = np.trapz(r * Pt,r )
    return Q
#########################################################################################
#########################################################################################
def getThrustFromBlade(r0, Pn0, R):
    n = len(r0)
    r = np.zeros((n + 1))
    Pn = np.zeros((n + 1))
    r[0:n] = r0
    r[-1] = R
    Pn[0:n] = Pn0
    T = np.trapz(Pn,r)
    return T
#########################################################################################
#########################################################################################
def fInductionCoefficients(a_last, Vrel_norm, Un, Ut, V0_in3, V0_in4, nW_in4, omega, chord, F, Ftip, CnForAI, CtForTI,
                           lambda_r, sigma, phi, Algo):
    a = 1. / ((4. * F * sin(phi * pi / 180.) ** 2) / (sigma * CnForAI) + 1.)
    # Thrust coefficient from the momentum theory => alast
    # CT=(1-a_last).^2.*sigma.*CnForAI./((sind(phi)).^2)
    # print Vrel_norm, sigma, CnForAI, np.linalg.norm(V0_in3)
    CT_loc = np.linalg.norm(Vrel_norm) ** 2 * sigma * CnForAI / (np.linalg.norm(V0_in3) ** 2)  # that's a CT loc
    # --------------------------------------------------------------------------------
    # --- Hight thrust correction
    # --------------------------------------------------------------------------------
    a, CT_loc = fCorrectionHighThrust(Algo.CTcorrection, a, CnForAI, phi, a_last, sigma, F, Ftip, CT_loc)
    a = a * Algo.relaxation + (1. - Algo.relaxation) * a_last

    if Algo.bSwirl is True:
        aprime = 1. / ((4. * F * sin(phi * pi / 180.) * cos(phi * pi / 180.)) / (sigma * CtForTI) - 1.)
        if Algo.SwirlMethod == 'HAWC':
            aprime = (np.linalg.norm(Vrel_norm) ** 2 * CtForTI * sigma) / (
                4. * (1. - a) * np.linalg.norm(V0_in4) ** 2 * lambda_r)
        else:
            raise Exception('Method not implemented')
    else:
        aprime = a * 0.

    return a, aprime, CT_loc
#########################################################################################
#########################################################################################
def fCorrectionHighThrust(CTcorrection, a, CnForAI, phi, a_last, sigma, F, Ftip, CT_loc):
    # print type(a)
    # print type(a) is not 'numpy.float64'
    # if type(a) == 'numpy.float64':
    #     raise Exception('Ctcorrection implemented for none array')
    if CTcorrection == 'GlauertCT':
        ac = 0.3
        if (a > ac):
            fg = 0.25 * (5. - 3. * a)
            a = CT_loc / (4. * F * (1. - fg * a))
    elif CTcorrection == False:
        pass
    else:
        raise Exception('CT correction not implemented')
    return a, CT_loc
#########################################################################################
#########################################################################################
def fAeroCoeffWrap(Rotor, PcDataAll, e, alpha, phi, chord, Vrel_norm, Re, Fperf, WT, Algo, Misc):
    fs = 0;    ne = 1
    Cl = [];    Cd = []

    if Misc.Format == 'flex':
        ee = int(Rotor.ProfileSet[1, e])
        # print ee
        Cd = np.interp(alpha, PcDataAll[ee - 1][:, 0], PcDataAll[ee - 1][:, 2])
        if Algo.bDynaStall == True:
            raise Exception('Not implemented')
        else:
            Cl = np.interp(alpha, PcDataAll[ee - 1][:, 0], PcDataAll[ee - 1][:, 1])
    else:
        ClCdCm = fAeroCoeff(alpha, PcDataAll, Rotor.ProfileSet[:, e], Rotor.Profile_thickness_rel,
                            Rotor.thickness_rel[e], Re, Algo.bReInterp, Algo.bThicknessInterp, Algo.bRoughProfiles)
        Cl = ClCdCm[0]
        Cd = ClCdCm[1]
        # print ClCdCm
    # Normal and tangential
    CnNoDrag = Cl * cos(phi * pi / 180.)
    CtNoDrag = Cl * sin(phi * pi / 180.)
    Cn = Cl * cos(phi * pi / 180.) + Cd * sin(phi * pi / 180.)
    Ct = Cl * sin(phi * pi / 180.) - Cd * cos(phi * pi / 180.)
    # performance correction on Cn Ct
    Cn = Fperf * Cn
    Ct = Fperf * Ct
    if (Algo.bAIDrag):
        CnForAI = Cn
    else:
        CnForAI = Fperf * CnNoDrag

    if (Algo.bTIDrag):
        CtForTI = Ct
    else:
        CtForTI = CtNoDrag
    return Cl, Cd, Cn, Ct, CnForAI, CtForTI, fs

def fAeroCoeff(alpha, PcDataAll, ProfileSet, Profile_thickness_rel, rel_thickness, Re, bReInterp, bThicknessInterp,
               bRough):
    if bRough:
        raise Exception('Rough profile are not implemented yet')
    else:
        Data = PcDataAll
    ClCdCm = np.zeros((3))
    temp = np.zeros((3, 3))
    ii1 = int(ProfileSet[1])
    ii2 = int(ProfileSet[2])
    if bReInterp:
        raise Exception('Reynolds interpolation not implemented yet')
    else:
        if ProfileSet[0] == 1:
            for j in np.arange(1, 4, 1):
                ClCdCm[j - 1] = np.interp(alpha, PcDataAll[ii1 - 1][:, 0], PcDataAll[ii1 - 1][:, j])
        else:
            # first profile
            for j in np.arange(1, 4, 1):
                temp[2, j - 1] = np.interp(alpha, PcDataAll[ii1 - 1][:, 0], PcDataAll[ii1 - 1][:, j])
            if bThicknessInterp is False:
                ClCdCm = temp[2, 0:2]
            else:
                # second profile
                for j in np.arange(1, 4, 1):
                    temp[1, j - 1] = np.interp(alpha, PcDataAll[ii2 - 1][:, 0], PcDataAll[ii2 - 1][:, j])
                for j in np.arange(1, 4, 1):
                    ClCdCm[j - 1] = np.interp(rel_thickness, Profile_thickness_rel[ii1 - 1:ii2], temp[2:0:-1, j - 1])
    return ClCdCm
#########################################################################################
#########################################################################################
def fSetRotorGrid(ngrid, Rotor, Env, WT, Spec, State, Algo, Misc):
    if (ngrid == 0) == 1:
        print 'You didn''t specify a ngrid parameter'
    else:
        rfull = np.linspace(Rotor.rhub, Rotor.R, ngrid + 1)
        r_mid = (rfull[0:-1] + rfull[1:]) / 2.
        Rotor = fInterpRotor(r_mid, Rotor)
    return Rotor
#########################################################################################
#########################################################################################
def fInterpRotor(r_mid, Rotor):
    Rotor.chord = np.interp(r_mid, Rotor.r, Rotor.chord)
    Rotor.thickness_rel = np.interp(r_mid, Rotor.r, Rotor.thickness_rel_prof)
    Rotor.twist = np.interp(r_mid, Rotor.r, Rotor.twist)
    Rotor.r = r_mid

    Rotor.ProfileSet = np.zeros((3, len(r_mid)))

    for i in np.arange(len(r_mid)):
        profilesup = np.where(Rotor.Profile_thickness_rel >= Rotor.thickness_rel[i])
        profileinf = np.max(np.where(Rotor.Profile_thickness_rel <= Rotor.thickness_rel[i]))

        if not profilesup:
            profilesup = profileinf
        elif not profileinf:
            profileinf = 0.
        profilesup = profilesup[0][0]
        cp = int(profileinf != profilesup)
        # b=np.array([cp+1, profileinf, profilesup])
        Rotor.ProfileSet[:, i] = np.array([cp + 1, profileinf + 1, profilesup + 1])

    Rotor.e_ref_for_khi = np.min(np.where(Rotor.r > 0.7 * Rotor.R))
    Rotor.ne = len(Rotor.r)
    ## Rotor dr if needed
    Rotor.dr = Rotor.r * 0.
    Rotor.dr[0] = 2.0 * (Rotor.r[0] - Rotor.rhub)
    for i in np.arange(1, len(Rotor.r)):
        Rotor.dr[i] = 2.0 * (Rotor.r[i] - Rotor.r[i - 1] - Rotor.dr[i - 1] * 0.5)
    return Rotor
#########################################################################################
#########################################################################################
def fInitWT(sWT, Format, pathtowt):
    Rotor = InitRotor()
    Env = InitEnv()
    WT = InitTurbine()
    Spec = InitSpec(Env, Rotor)
    #    Controller=InitController()
    State = InitState()
    Algo = InitAlgo()
    Misc = InitMisc()
    parMisc = {'WTname': sWT, 'Format': Format}
    Misc.Set(**parMisc)

    HtcFile = glob(pathtowt + sWT + '/*.htc')
    #    AeFile= glob(pathtowt+sWT+'/data/*ae*')
    #    PcFile= glob(pathtowt+sWT+'/data/*pc*')
    SpecFile = glob(pathtowt + sWT + '/*Spec*')
    Rotor, Spec, WT, Algo = fReadSpec(SpecFile, WT, Rotor, Spec, Algo)
    AeSet, PcFileName, AeFileName, Nb, PitchAxis = fReadHtcFile(HtcFile, 'blade1')
    AeData = fReadAeFile(AeFileName, pathtowt, sWT, AeSet, 4)
    PcSet = AeData[0, 3]
    PcDataAll, thickness_rel_prof, ndata = fReadPcFile(PcFileName, pathtowt, sWT, PcSet)

    Rotor.r = AeData[:, 0] + Rotor.rhub
    Rotor.chord = AeData[:, 1]
    Rotor.thickness_rel_prof = AeData[:, 2]
    Rotor.Profile_thickness_rel = thickness_rel_prof
    if Format == 'hawc':
        Stations = PitchAxis
        twist = Stations[:, 3]
        rtwist = Stations[:, 2] + Rotor.rhub
    # Dealing with sign
    if (np.mean(twist) < 0) == 1:
        sign = -1;
    else:
        sign = 1
    # Dealing with problem of interpolation
    if (max(rtwist) < max(Rotor.r)) == 1:
        print '! For twist interpolation, last blade section in htc file should be at bigger the last of the ae file. Im replacing it....'
        rtwist[-1] = max(Rotor.r)

    if (min(rtwist) > min(Rotor.r)) == 1:
        print '! For twist interpolation, first blade section in htc file should be at smaller (usually 0) than the one in the ae file. Im replacing it....'
        rtwist[0] = Rotor.r[0]

    # Interpolating twist
    Rotor.twist = np.interp(Rotor.r, rtwist, twist) * sign

    Rotor.R_coned = Rotor.R * cos(Rotor.cone * pi / 180.)  # Rotor coned radius [m]
    Rotor.SweptArea = pi * (Rotor.R * cos(Rotor.cone * pi / 180.)) ** 2
    Spec.TSR_rated = Spec.Omega_rated * Rotor.R * cos(Rotor.cone * pi / 180.) / Spec.V_rated
    return PcDataAll, AeData, Rotor, Env, WT, Spec, State, Algo, Misc
#########################################################################################
#########################################################################################
def fReadPcFile(PcFileName, pathtowt, sWT, PcSet):
    fd = open(pathtowt + sWT + PcFileName[1:], 'r')
    line = fd.readline()
    n = line.split(" ");
    n = filter(None, n)
    #    NrSet = float(n[0])
    for i in np.arange(0, PcSet - 1, 1):
        line = fd.readline()
    line = fd.readline()
    n = line.split(" ");
    n = filter(None, n)
    NrSubSet = int(n[0])
    PcData = np.zeros((NrSubSet, 4))
    PcDataAll = [i for i in range(NrSubSet)]
    thickness_rel_prof = np.zeros((NrSubSet))
    ndata = np.zeros((NrSubSet))
    for i in np.arange(0, NrSubSet, 1):
        line = fd.readline()
        n = line.split(" ");
        n = filter(None, n)
        Tempvec = n[2:0:-1]
        Tempvec = [float(l) for l in Tempvec]
        PcData = np.zeros((len(np.arange(0, Tempvec[1], 1)), 4))
        thickness_rel_prof[i] = Tempvec[0]
        ndata[i] = Tempvec[1]
        for j in np.arange(0, Tempvec[1], 1):
            line = fd.readline()
            n = line.split(" ");
            n = filter(None, n)
            PcData[j, :] = n[0:4]
        PcDataAll[i] = PcData
    fd.close()
    return PcDataAll, thickness_rel_prof, ndata
#########################################################################################
#########################################################################################
def fReadAeFile(AeFileName, pathtowt, sWT, AeSet, ncol):
    fd = open(pathtowt + sWT + AeFileName[1:], 'r')
    line = fd.readline()
    n = line.split(" ");
    n = filter(None, n)
    NrSet = float(n[0])
    for i in np.arange(0, AeSet[0] - 1, 1):
        line = fd.readline()
    line = fd.readline()
    n = line.split(" ");
    n = filter(None, n)
    Label = n[-1]
    Nr = float(n[1])
    AeData = np.zeros((Nr, ncol))
    if len(Label) == 0:
        Label = None
    for i in np.arange(0, Nr, 1):
        line = fd.readline()
        n = line.split("\t"); n = [w.replace('\n', '') for w in n]
        n = filter(None, n)
        AeData[i, :] = [float(k) for k in n]
    return AeData
#########################################################################################
#########################################################################################
def fReadHtcFile(HtcFile, BladeBodyName, ):
    fd = open(HtcFile[0], 'r')
    while 1:
        line = fd.readline()
        if not line:
            break
        pass
        ## READ PITCH AXIS DATA
        if 'name' in line and BladeBodyName in line:
            while 1:
                line = fd.readline()
                if 'begin' in line and 'c2_def' in line:
                    break
            line = fd.readline()
            line=line.replace('\t', ' ')
            n = line.split(" ");
            n = filter(None, n)
            nsec = float(n[n.index('nsec') + 1])
            PitchAxis = np.zeros((nsec, 4))
            for i in np.arange(0, float(nsec), 1):
                line = fd.readline()
                line=line.replace('\t', ' ')
                n = line.split(' ')
                n = filter(None, n)
                nn = [float(j) for j in n[2:6]]
                PitchAxis[i] = nn
                ## READ AERODYNAMIC FILE
        if 'begin' in line and 'aero ' in line:
            while not 'end aero' in line:
                line = fd.readline()
                if 'nblades' in line:
                    n = line.split(" ");
                    n = filter(None, n)
                    n = [w.replace('\t', '') for w in n]
                    Nb = float(n[n.index('nblades') + 1][0])
                if 'ae_filename' in line:
                    n = line.split(" ");
                    n = filter(None, n)
                    n = [w.replace('\t', '') for w in n]
                    AeFileName = n[n.index('ae_filename') + 1]
                if 'pc_filename' in line:
                    n = line.split(" ");
                    n = filter(None, n)
                    n = [w.replace('\t', '') for w in n]
                    PcFileName = n[n.index('pc_filename') + 1]
                if 'ae_sets' in line:
                    n = line.split(" ");
                    n = filter(None, n)
                    n = [w.replace('\t', '') for w in n]
                    AeSet = n[n.index('ae_sets') + 1:-1]
                    AeSet = [float(i) for i in AeSet]
    fd.close()
    return AeSet, PcFileName, AeFileName, Nb, PitchAxis
#########################################################################################
#########################################################################################
def fReadSpec(SpecFile, WT, Rotor, Spec, Algo):
    fd = open(SpecFile[0], 'r')
    WsRpmPitch=[]
    A=[]
    while 1:
        line = fd.readline()
        if not line:
            break
        pass
        if 'SPECS' in line:
            # Read RPM and Pitch curve
            while 1:
                line = fd.readline()
                if not line:
                    break
                else:
                    n = line.split(" ");
                    n = filter(None, n)
                    nn = [float(j) for j in n[0:3]]
                    WsRpmPitch.append(nn)
        else: # read other specs
            n = line.split(" ");
            n = filter(None, n)
            A.append(n[0])
    fd.close()
    WsRpmPitch=np.array(WsRpmPitch).T
    A=np.array(A)
    A = [float(j) for j in A]
    parRotor = {'nB': A[0],
                'BladeLength':A[1],
                'rhub': A[2],
                'cone':A[3],
                'R': A[1] + A[2]}
    Rotor.Set(**parRotor)

    parWT = {'tilt': A[4],
             'yaw': A[5],
             'H': A[6]}
    WT.Set(**parWT)

    Algo.Ngrid = A[7]
    parAlgo = {'Ngrid': A[7]}
    Algo.Set(**parAlgo)

    parSpec = {'Omega_rated': A[8] * 2. * pi / 60.,
               'P_rated': A[9],
               'Pitch': WsRpmPitch[2,:],
               'RPM': WsRpmPitch[1,:],
               'WS': WsRpmPitch[0,:]}
    Spec.Set(**parSpec)

    parRotor = {'Omega': Spec.Omega_rated}
    Rotor.Set(**parRotor)
    return Rotor, Spec, WT, Algo
#########################################################################################
#########################################################################################
