# -*- coding: utf-8 -*-
"""
Created on Tue May  5 08:44:05 2015

@author: Ewan Machefaux
# Python interpreation of Matlab-based steady Blade Element Momentum of Emmanuel Branlard (emmanuel.branlard@gmail.com)
"""

# import matplotlib as plt
import time
from sys import exit
import numpy as np
from math import cos, atan2, acos, exp, sqrt, sin, pi
from cBEM import *
from os import chdir, listdir, path, getcwd
from glob import glob
from DWM_misc import find_intersections, unique_rows, unique2d
import os
# import matlab.engine
import matplotlib.pyplot as plt

from scipy.optimize import fmin_cobyla

#from Read_Flex5 import read_Pitch_Flex5


#########################################################################################
#########################################################################################
def getInduction(ngrid, sWT, Format, U0, meta, ID_waked, **kwargs):
    """Function that compute the axial induction distribution for subsequent use in the Ainslie wake flow model

        Inputs
        ----------
        ngrid (int): number of grid point per radius in BEM
        sWT (string): turbine name
        Format (string) : input file format (HAWC2)
        U0: mean hub height velocity
        meta (instance of class): Instance of class Meta holding DWM core variables

        Outputs
        ----------
        BEM (instnce of class) : refer to class description for details
    """

    opt = 'python'
    # opt = 'matlab'
    derating = kwargs.get('derating')
    # Initialization
    PcDataAll, AeData, Rotor, Env, WT, Spec, State, Algo, Misc = fInitWT(sWT, Format, '../WT-data/')
    Rotor = fSetRotorGrid(ngrid, Rotor, Env, WT, Spec, State, Algo, Misc)
    # if derating == 1.:

    RPM=np.interp(U0,Spec.WS,Spec.RPM)
    Pitch=np.interp(U0,Spec.WS,Spec.Pitch)         # PITCH control from basic sdwm data

    #For Comparison with Flex5
    """
    if ID_waked == {'1': [], '0': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': []}:
        print 'Calculation for 1st Turbine in the row is based on Flex5 data for Pitch and RPM'
        WTG = 'NREL5MW'
        INPUT = 'sweep'
        dir = 'C:/Users/augus/OneDrive/Documents/Stage/Results/Flex5/' + WTG + '/' + INPUT + '/Data/'
        kind = 'tim'

        DATA = np.loadtxt(dir + 'converged_' + kind + '.dat')  # (s), (m/s), (deg), (kW), (...)...

        RPM=DATA[3,3] * 30./pi                             # Common RPM between sDWM BEM and BEM
        Pitch = np.interp(U0, DATA[:,1], DATA[:,2])          # PITCH control from Flex5 data for comparison with flex5
    #"""


    # RPM, Pitch= getOptimalLambdaPitch(meta,U0,derating,Spec,Rotor)
    # float("{0:.2f}".format(meta.WTG_spec.get_P(meta.mean_WS_DWM)))


    Sim = InitSim()
    parSim = {'Name': sWT,
          'WT': sWT,
          'WS': U0,
          'RPM': RPM,
          'PITCH': Pitch,
          #'RPM_opt' : np.interp(U0,Spec.WS,Spec.RPM),
          #'PITCH_opt' : np.interp(U0,Spec.WS,Spec.Pitch)
          'RPM_opt': RPM,
          'PITCH_opt' : Pitch
            }
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
    #print 'float("{0:.2f}".format(meta.WTG_spec.P_rated *derating))',float("{0:.2f}".format(meta.WTG_spec.P_rated *derating))
    # print 'float("{0:.2f}".format(BEM.Power/1000.)', float("{0:.2f}".format(BEM.Power/1000.))
    ################################################################
    ################################################################
    return BEM
#########################################################################################
#########################################################################################


def fBEMsteady(WT, Sim, Wind, Algo, Rotor, PcDataAll, Env, Spec, State, Misc):
    """Function that run the BEM main loop

        Inputs
        ----------
        Sim (instance of class) holds simulation related parameters
        Wind (instance of class) holds wind inflow related parameters
        Algo (instance of class) holds BEM algorithm related parameters
        Rotor (instance of class) holds rotor related geometrical parameters
        PcDataAll (array of float) holds the profile coefficient Cm, Cl, Cd and AoA
        Env (instance of class) holds ambient conditions related parameters
        Spec (instance of class) holds turbine specification parameters
        Misc (instance of class) holds I/O parameters

        Outputs
        ----------
        BEM (instnce of class) : refer to class description for details
    """
    cone = Rotor.cone
    # print 'cone is', Rotor.cone
    V0 = Wind.V0[2] * cos(cone * pi / 180.)
    # V0 = Wind.V0[2]
    VHubHeight = Wind.V0[2]
    nB = Rotor.nB
    ne = Rotor.ne
    r = Rotor.r
    #print r
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
                if Algo.TipLossMethod == 'Prandtl':     # originally written Glauert but it's Prandtl
                    #print "Performing Prandtl's Tip Loss Factor"
                    if (sin(phi * pi / 180.) > 0.01):
                        # prandtl tip correction
                        Ftip = 2. / pi * acos(exp(-nB / 2. * (R - r[e]) / (r[e] * sin(phi * pi / 180.))))
                #elif Algo.TipLossMethod == 'Prandtl':    # I think it's not Prandlt
                else:     #Where does it come from? Is it another formula for Prandtl? I don't think so (no relation with phi)
                    Ftip = 2. / pi * acos(exp(-nB / 2. * (1. - lambda_r[e] / lambda_) * sqrt(1. + lambda_ ** 2)))
            F = Ftip
            # here will be implemented Hub losses in the future
                    # Implemented by augr
            if Algo.bHubLoss:
                # q = B / 2. * (r[k]-R_end_cylinder_according_chord) / (R_end_cylinder_according_chord * np.sin(phi[k])) # my BEM
                Fhub = 2. / pi * acos(exp(-nB / 2. * (r[e] - rhub) / (rhub * sin(phi * pi / 180.))))
                if Algo.bTipLoss:
                    F = Fhub * Ftip
                else:
                    F = Fhub
            # --------------------------------------------------------------------------------
            # --- Step 3: Angle of attack
            # --------------------------------------------------------------------------------
            #pitch=pi/2             #lock the pitch angle, otherwise put in comment
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
    #print 'BEM.RPM', BEM.RPM

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
    """Function that estimates the aerodynamic torque

        Inputs
        ----------
        r0: discrete radial position of blade [m]
        Pt0: tangential load in N
        R: blade radius [m]

        Outputs
        ----------
        Q: aerodynamic torque [Nm]
    """
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
    """Function that estimates the thrust force

        Inputs
        ----------
        r0: discrete radial position of blade [m]
        Pn0: normal load in N
        R: blade radius [m]

        Outputs
        ----------
        T: aerodynamic thrust [N]
    """
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
    """Function that compute the induction coefficients while applying swirl

        Inputs
        ----------
        a_last (float) : last iteration axial induction factor
        Vrel_norm (float): normed relative velocity
        Un (float): normal velocity
        Ut (float): tangential velocity
        V0_in3 (float): [0., 0., V0]
        V0_in4 (float): [0., 0., V0]
        nW_in4 (float): [0., 0., -a_last * V0]
        omega (float): rotor rotational speed rad/s
        chord (float): rotor chord length distribution
        F (float): total loss
        Ftip (float): tip loss
        CnForAI (float): normal force coefficient
        CtForTI (float): tangential force coefficient
        lambda_r (float): speed ratio distribution
        sigma (float): blade solidity
        phi (float): flow angle
        Algo (instance of class): holds algorith BEM related parameters

        Outputs
        ----------
        a: axial induction factor
        aprime: tangential induction factor
        CT_loc: local thrust coefficient
    """
    a = 1. / ((4. * F * sin(phi * pi / 180.) ** 2) / (sigma * CnForAI) + 1.)
    # Thrust coefficient from the momentum theory => alast
    # CT=(1-a_last).^2.*sigma.*CnForAI./((sind(phi)).^2)
    #print Vrel_norm, sigma, CnForAI, np.linalg.norm(V0_in3)
    CT_loc = np.linalg.norm(Vrel_norm) ** 2 * sigma * CnForAI / (np.linalg.norm(V0_in3) ** 2)  # that's a CT loc
    #print 'Ct_loc', CT_loc
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
    """Function that applying high thrust coefficient correction

        Inputs
        ----------
        CTcorrection (string) : type of high thrust coefficient
        a (float): axial induction coefficient
        CnForAI: normal force coefficient
        phi: flow angle
        a_last: last iteration axial induction
        sigma: blade solidity
        F: total losses
        Ftip: tip losses
        CT_loc: local thrust coefficient

        Outputs
        ----------
        a: axial induction factor
        CT_loc: local thrust coefficient
    """
    # print type(a)
    # print type(a) is not 'numpy.float64'
    # if type(a) == 'numpy.float64':
    #     raise Exception('Ctcorrection implemented for none array')
    if CTcorrection == 'GlauertCT':
        #print 'Performing Glauert Correction for Ct'
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
    """Tabulated airfoil data interpolation

        Inputs
        ----------
        Rotor (instance of class) holds rotor related parameter
        PcDataAll (array of float) holds the profile coefficient Cm, Cl, Cd and AoA
        e : element index in BEM loop
        alpha: Angle Of Attack
        phi: flow angle
        chord: chord length
        Vrel_norm: normed relative velocity
        Re: Reynolds number
        Fperf: performance losses (not implemented)
        Algo (instance of class) holds alogrithm BEM related parameter
        Misc (instance of class) holds I/O parameters

        Outputs
        ----------
        Cl: lift coefficient
        Cd: drag coefficient
        Cn: normal coefficient
        Ct: tangential coefficient
        CnForAI: normal coefficient for correction
        CtForTI: tangential coefficient for correction
    """
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
    """Function that interpolates in the profile coefficients tabulated airfoil data

        Inputs
        ----------
        alpha: angle of attack deg
        PcDataAll: array holding all profile coefficients
        ProfileSet: index of profile set coefficients for interpolation
        Profile_thickness_rel: relative thickness of each of the profile coeff sets
        rel_thickness: rel thickness at the point of interpolation
        Re: Reynolds number
        bReInterp: Reynolds number interpolation flag
        bThicknessInterp: Thickness interpolation flag
        bRoughs: Roughness profile interpolation flag

        Outputs
        ----------
        ClCdCm (float): vector containing the Cl, Cd and Cm values at interpolation point
    """

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
    """Function that discretizes the blade in elements

        Inputs
        ----------
        ngrid (int): numnber of blade elements-1
        Rotor (instance of class) holds rotor related geometrical parameters
        Outputs
        ----------
        Rotor (instance of class): updated instance of class Rotor
    """
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
    """Function that interpolate the blade geometry at the blade element center

        Inputs
        ----------
        r_mid float): blade element center
        Rotor (instance of class) holds rotor related geometrical parameters
        Outputs
        ----------
    Rotor (instance of class): updated instance of class Rotor
    """
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
    """Function that initializes the input reader
    """

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
    #print 'Rotor.cone: ' ,Rotor.cone
    Rotor.SweptArea = pi * (Rotor.R * cos(Rotor.cone * pi / 180.)) ** 2
    Spec.TSR_rated = Spec.Omega_rated * Rotor.R * cos(Rotor.cone * pi / 180.) / Spec.V_rated
    return PcDataAll, AeData, Rotor, Env, WT, Spec, State, Algo, Misc
#########################################################################################
#########################################################################################
def fReadPcFile(PcFileName, pathtowt, sWT, PcSet):
    """Function that loads HAWC2 profile coefficients file
    """
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
        #print 'Tempvec float: ', Tempvec
        Tempvec = [int(l) for l in Tempvec]
        #print 'Tempvec int: ', Tempvec
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
    """Function that loads HAWC2 aerodynamic coefficients file
    """
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
    #print 'ncol: ', ncol
    #print 'Nr: ', Nr
    Nr = int(Nr)
    #print 'Nr: ', Nr
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
    """Function that loads HAWC2 HTC file
    """
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
            #print 'nsec float: ', nsec
            nsec = int(n[n.index('nsec') + 1])
            #print 'nsec int: ', nsec
            PitchAxis = np.zeros((nsec, 4))
            for i in np.arange(0, nsec, 1):
                line = fd.readline()
                line=line.replace('\t', ' ')
                n = line.split(' ')
                n = filter(None, n)
                nn = [float(j) for j in n[2:6]]
                #print 'nn float: ', nn
                #nn = [int(float(j)) for j in n[2:6]]
                #print 'nn int: ', nn
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
    """Function that loads the turbine spec file
    """
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
    """
    PITCH = WsRpmPitch[2,:]
    RPM = WsRpmPitch[1,:]
    WS = WsRpmPitch[0,:]

    plt.figure(), plt.title('Wind Turbine Data'), plt.plot(WS,PITCH,label='Pitch(U)'), plt.legend(),plt.show()
    plt.figure(), plt.title('Wind Turbine Data'), plt.plot(WS, RPM, label='Omega(U)'), plt.legend(), plt.show()
    #"""
    parRotor = {'Omega': Spec.Omega_rated}
    Rotor.Set(**parRotor)
    return Rotor, Spec, WT, Algo
#########################################################################################
#########################################################################################