# -*- coding: utf-8 -*-
""" Multiple classes definition from the Blade Element Momentum
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
This class structure follows the Matlab implementation of Dr Emmnanuel Branlard
<emmanuel.branlard@gmail.com>. The present steady BEM implementation is a very limited part
of the Matlab BEM, therefore some classes are limited in size as only steady state related variables
are used. Future development will include unsteady BEM implementation
with non symmetric inflow conditions (shear and veer) as well as blade structural properties.
"""

########################################################################
### Rotor definition
########################################################################
## Main parameters 
from math import gamma, pi
import numpy as np


class InitRotor:
    """ This class holds rotor related parameters
    """
    def __init__(self):
        # Required Params
        self.nB = -1.0  # number of blades
        self.R = -1.0  # Rotor radius [m]
        self.BladeLength = -1.0  # Rotor radius [m]
        self.rhub = 0.0  # Hub  length [m]
        self.cone = 0.0  # cone angle [deg]
        self.ne = -1.0
        # self.M=-1;      # mass [kg]
        # self.I=-1;       #Inertia [kg/m/m]??? kg m^2
        self.rb_center_in4 = [0., 0., 0.]
        self.HubHeight = 60.0
        self.Omega = -1.  # Nominal Rotational speed [rad/s]

    def Set(self, **parRotor):
        """Parsing data ."""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)


class InitEnv:
    """ This class holds ambient conditions parameters
    """
    def __init__(self):
        self.rho = 1.225  # air density [kg/m^3]
        self.g = 9.81  # gravity [m/s^2]
        self.KinVisc = 15.68 * 10 ** -6  # Knematic viscosity nu [m^2/s]
        self.A = 6.67  # shape factor # from 31784
        self.k = 2.04  # scale factor

    def Set(self, **parEnv):
        """Parsing data ."""
        for key, value in parEnv.iteritems():
            setattr(self, key, value)


class InitTurbine:
    def __init__(self):
        # Required Params
        ## Shaft
        self.Lshaft = 0.0  # Shaft length [m]
        self.rs_in2 = [0.0, 0.0, -self.Lshaft]
        # Shaft.k=-1; # [Nm/rad] 4*10^7
        ## Generator
        self.I = -1.0  # 7.95*10^5 [kg/m^2]
        self.fMoment = 2 * 10 ** 6  # 2*10^6 *(max(x,Rotor.Omega)-Rotor.Omega);
        ## Nacelle
        # Nacelle.M=-1;  # [kg]
        self.tilt = 0  # tilt angle [deg]
        ## Tower
        self.H = 60.0  # height of the tower [m]
        self.r1 = 2.125
        self.r2 = 2.375
        self.H1 = self.H / 2.0
        self.Xtower = [-self.r2, -self.r2, -self.r1, self.r1, self.r2, self.r2]
        self.Ytower = [0.0, self.H1, self.H, self.H, self.H1, 0.0]
        self.rt_in1 = [self.H, 0.0, 0.0]
        # Tower.k=-1; # [N/rad] 2*10^6

    def Set(self, **parTurbine):
        """Parsing data ."""
        for key, value in parTurbine.iteritems():
            setattr(self, key, value)


class InitSpec:
    """ This class holds turbine operating specifications
    """
    def __init__(self, Env, Rotor):
        # Required Params
        self.Cp_rated = -1.0  # estimation of Cp_max at mean wind speed from 31784
        self.P_rated = -1.0  # 500kW
        self.V_rated_thumbs = gamma(1. + 1. / Env.k) * Env.A + 6.  # rule of thumb mean WS + 6m/s
        self.V_rated = (self.P_rated / (self.Cp_rated * 0.5 * Env.rho * Rotor.R ** 2. * pi)) ** (1. / 3.)
        self.Reynolds_inf = 0.  # 75000*Rotor.R
        self.Reynolds_sup = 0.  # 150000*Rotor.R
        self.Omega_rated = 0.
        self.TSR_rated = 0. # tip speed ratio at rated
        self.RPM = [] # rotational speed in RPM
        self.Pitch = [] # pitch angle in deg
        self.WS = [] # mean hub height wind speed
    def Set(self, **parSpec):
        """Parsing data ."""
        for key, value in parSpec.iteritems():
            setattr(self, key, value)


class InitController:
    def __init__(self):
        # Required Params
        self.fpitch = 0.  #
        self.fyaw = 0.  # 
        self.lastpitch = 0.  #

    def Set(self, **parController):
        """Parsing data ."""
        for key, value in parController.iteritems():
            setattr(self, key, value)


class InitState:
    def __init__(self):
        # Required Params
        self.pitch = 0.  # [deg]
        self.yaw = 0.  # [deg]
        self.psi = 0.  # [deg]
        self.w_guess = -2.
        #### HERE WE SHOULD MAKE A SUBCLASS LAST !!
        #        self.last.W=0.     #ones(3,ne,nB).*w_guess;       #temporary induced velocity
        #        self.last.W0=0.    #ones(3,ne,nB).*w_guess;       #temporary induced velocity
        #        self.last.W_qs=0.  #ones(3,ne,nB).*w_guess; #temporary quasistatic induced velocity
        #        self.last.W_int=0. # ones(3,ne,nB).*w_guess; #temporary intermediate induced velocity
        #        self.last.fs=0.
        self.Power = 0.
        self.lastPower = 0.

    def Set(self, **parState):
        """Parsing data ."""
        for key, value in parState.iteritems():
            setattr(self, key, value)


class InitMisc:
    """ This class holds turbine operating initial states
    """
    def __init__(self):
        # Required Params
        self.Profiles = -1. # profile sets  in AE HAWC2 file
        self.Format = 'hawc' # format of turbine files
        self.WTname = 'NREL5MW' # turbine name for input files

    def Set(self, **parMisc):
        """Parsing data ."""
        for key, value in parMisc.iteritems():
            setattr(self, key, value)


class InitAlgo:
    """ This class holds BEM algorithm flags and parameters
    """
    def __init__(self):
        self.nbIt = 200  # maximum number of iterations in BEM
        self.aTol = 10 ** -6 # tolerance for axial induction factor convergence
        self.relaxation = 0.5  # relaxation factor in axial induction factor
        self.CTcorrection = 'GlauertCT'  #  type of CT correction more model implementated in the future like 'spera'
        self.Ngrid = 1.0
        self.bSwirl = True  # swirl flow model enabled / disabled
        self.SwirlMethod = 'HAWC' # type of swirl model
        self.bTipLoss = True # enable / disable tip loss model
        self.bHubLoss = True # enable / disable hub loss model
        self.bTipLossCl = False # enable / disable Cl loss model
        self.TipLossMethod = 'Prandtl'  # type of tip loss model  # I change it to Prandtl (originally Glauert)
        self.bDynaStall = 0 # dynamic stall model
        self.bAIDrag = True # influence on drag coefficient on normal force coefficient
        self.bTIDrag = True # influence on drag coefficient on tangential force coefficient
        self.bReInterp = False # interpolate the input tabulated airfoil data for Reynolds variation
        self.bThicknessInterp = True # interpolate the input tabulated airfoil data for thickness variation
        self.bRoughProfiles = False # use rough profiles for input airfoil data

    def Set(self, **parAlgo):
        """Parsing data ."""
        for key, value in parAlgo.iteritems():
            setattr(self, key, value)


class InitSim:
    """ This class holds further simulation parameters
    """
    def __init__(self):
        self.WT = '' # turbine name
        self.Name = '' # turbine name, redundant to Name, adaptation from Matlab library
        self.rho = 1.225 # air density
        self.KinVisc = 15.68 * 10 ** -6 # kinematic viscosity
        self.WS = 0. # wind speed m/s
        self.RPM = 0. # rotational speed RPM
        self.PITCH = 0. # pitch angle [deg]
        self.YAW = 0. # yaw angle [deg]
        self.Omega = 0. # rotational speed rad/s

    def Set(self, **parSim):
        """Parsing data ."""
        for key, value in parSim.iteritems():
            setattr(self, key, value)


class InitWind:
    """ This class holds inflow wind model. Further implementation will include non symmetric inflow
    (shear and veer) as well as inflow turbulence
    """
    def __init__(self):
        self.nu = 0.0
        self.model = 'Constant'
        self.V0 = 0.

    def Set(self, **parWind):
        """Parsing data ."""
        for key, value in parWind.iteritems():
            setattr(self, key, value)


class InitBEM:
    """ This class holds outputs from the BEM simulations
    """
    def __init__(self, ne):
        """ Initializing the data ."""
        self.phi = np.zeros((ne)) # flow angle in deg
        self.alpha = np.zeros((ne)) # angle of attack in deg
        self.a = np.zeros((ne)) # axial induction
        self.a_last = np.zeros((ne)) # last iteration axial induction for iterative loop
        self.aprime = np.zeros((ne)) # tangential induction
        self.aprime_last = np.zeros((ne)) # last iteration tangential induction
        self.Cd = np.zeros((ne)) # drag coefficient
        self.Cl = np.zeros((ne)) # lift coefficient
        self.Cn = np.zeros((ne)) # normal force coefficient
        self.Ct = np.zeros((ne)) # tangential force coefficient
        self.CT = np.zeros((ne)) # thrust coefficient
        self.Vrel = np.zeros((ne)) # relative velocity in m/s
        self.Re = np.zeros((ne)) # Reynolds number based on chord
        self.F = np.zeros((ne)) # loss corrections
        self.Fperf = np.zeros((ne))
        self.Fshen = [] # Shen's model loss corrections
        self.Un = np.zeros((ne)) # normal velocity
        self.Ut = np.zeros((ne)) # tangential velocity
        self.nIt = np.zeros((ne)) # number of iterations
        self.omega = [] # rotational speed in rad/s
        self.derated=False
        self.Ud=[] # demanded equivalent rotor averaged velocity from 
    def Set(self, **parBEM):
        """Parsing data ."""
        for key, value in parBEM.iteritems():
            setattr(self, key, value)

# class InitAerodynamic:
#    def __init__(self):
#        self.V0=[]
#        self.nu=[]
#        self.Model=[]
#        self.last=[]
#    def Set(self,V0,nu,Model,last):
#        self.V0=[0, 0, V0] #[0, 0, Vz] # wind speed Vx, Vy, Vz [m/s]
#        self.nu=nu #0.2 # shear coeff
#        self.Model=Model #'Constant' # turbulent or constant 
#        self.last=last #0
