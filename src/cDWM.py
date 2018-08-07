# -*- coding: utf-8 -*-
""" Multiple classes definition from the DWM flowfield model
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""
import numpy as np
from math import pi

class Meta:
    """ Main class holding all meta data for sDWM core model
    """
    def __init__(self):
        """ Initializing the data / Set default value ."""
        self.previous_sDWM = True
        # ---------------------------- # PLOT SETTINGS OPTIONS # ----------------------------------------------------- #
        # Deficit Process
        self.BEM_plot = False
        self.BC_detail_plot = False
        self.Deficit_Process_Detail = False
        self.AINSLIE_EV_details = False

        # Meandering Process (more settings in cTurb)
        self.MEANDERING_plot = False
        self.MEANDERING_detail_plot = False
        self.MEANDERING_WS_plot = True
        self.MEANDERING_WS_added_plot = True
        self.MEANDERING_TI_plot = False

        # Wake Added Turbulence Process
        self.WaT_detail = False

        # Build-up process and result in term of global deficit and Turbulence
        self.TI_Dynamic_for_Loads_plot = False  # According to Keck Atmospheric shear and wake ... 2013-2015

        self.DEFICIT_plot = True
        self.DEFICIT_details = False # For Dynamic, it gives result in time


        # ------------------------------------------------------------------------------------------------------------ #
        ## Default Ambient conditions
        self.WS                    = 8.0 # wind speed m/s
        self.TI                    = 0.07 # turbulence intensity % X10 ?

        ## Default Park settings:
        self.WTG                   = 'V80'
        # self.aSource               = 'BEM'# 'Table' no longer available
        # self.WTG_R                 = 40
        # self.H_hub                 = 70.0

        ## Domain settings:
        self.lx                    = 3.2 # global domain width in R
        self.ly                    = 3.2 # global domain height in R
        self.hub_x                 = self.lx/2. # position of the rotor center
        self.hub_y                 = self.ly/2.
        self.hub_z                 = np.asarray([]) # Turbine spacing [D]
        self.dR                    = 40 # points per R in radial direction default 20
        self.dx                    = 20 # points per R in lateral direction  default 20
        self.dy                    = 20 # points per R in longitudinal direction  default 20
        self.dz                    = 10 # points per D in axial direction default 10
        self.nx                    = int(self.lx * self.dx)  # nb points in x (lateral)  direction global flow field
        self.ny                    = int(self.ly * self.dy)  # nb points in x (lateral)  direction global flow field
        self.nz                    = np.asarray([])  # nb points in z (streamwise) direction global flow field
        self.nt = 0
        self.time = []
        self.x_vec                 = np.arange(1,self.nx+1,1)*(1.0/self.dx) # grid coordinates in lat. direction [R]
        self.y_vec                 = np.arange(1,self.ny+1,1)*(1.0/self.dy) # grid coordinates in long. direction [R]
        self.z_vec                 = np.asarray([]) # grid coordinates in streamwise. direction [D]
        self.z_vec_old             = np.asarray([]) # grid coordinates in streamwise. direction [D] at previous turbine
        self.x_mat                 = np.tile(self.x_vec.reshape(len(self.x_vec),1), self.ny)  # matrix grid holding x coordinates
        self.y_mat                 = np.tile(self.y_vec,(self.nx,1)) # matrix grid holding y coordinates
        self.z_mat                 = np.asarray([]) # matrix grid holding z coordinates
        self.nr                    = round(np.sqrt((max(self.x_vec - self.hub_x))**2 + (max(self.y_vec - self.hub_y))**2 ) * self.dR) # number of position in radial coordinate system
        self.r_dist_2              = np.sqrt((self.x_mat - self.hub_x)**2 + (self.y_mat - self.hub_y)**2 )  #Find distance to centre of wake plane
        self.C2C                   = np.asarray([]) #the "Center to Center" distance between hubs

        # Eddy viscosity calc mixL model specific parameter
        self.R_WTG                 = 1.0 # normalized turbine radius for eddy-viscosity calculation
        self.D_WTG                 = 2.0 * self.R_WTG # normalized turbine radius for eddy-viscosity calculation
        self.U0                    = 1.0 # normalized free stream velocity for eddy-viscosiy calculation
        self.dz_mixl               = self.D_WTG / self.dz  #note  is per D for eddy-viscosity calculation
        self.dr_mixl               = self.dR # points per R in radial direction identical to global flow field
        self.dr                    = self.R_WTG / self.dr_mixl  #note is per R for eddy-viscosity calculation
        self.lr_mixl               = self.lx+1. # domain size in radial direction [R] (Should be about 0.5R wider flow field domain width)
        self.vr_mixl               = np.linspace(0,self.lr_mixl-self.dr,(self.lr_mixl)/self.R_WTG*self.dr_mixl)  # coordinate vector
        self.vr_m                  = np.arange(1,self.dr_mixl*self.lr_mixl+1,1,dtype=float) / self.dr_mixl  # coordinate vector for polar discretization
        self.dA_DWM                = np.concatenate(([0], (pi*self.vr_m[1:]**2 - pi*self.vr_m[0:-1]**2)), axis=0)  # area vector for polar discretization
        self.lz_mixl               = np.asarray([]) # in R: 1R longer than global flow field due to backward diff scheme, turbine position dependent

        #if self.Keck:
        ## Model inputs calibrated constant
        self.Model                 = 'mixL' # Eddy viscosity model as per Keck et al. [2]
        #self.fU                    = 1.10 # fU = 1.1 & fR = 0.98 yields the best fit to AL calibration data according to Keck et al. [5].
        #self.fR                    = 0.98 #
        self.fU                    = 1. #Madsen et al [2].
        self.fR                    = 1.
        self.atmo_stab             = 'N'
        self.k1                    = 0.0919 # k1 = 0.0919 & k2 = 0.0178 yields the best fit to AL calibration data according to Keck et al. [5].
        self.k2                    = 0.0178  # in Keck we can compare with k1 0.914 and k2 = 0.0216

        ################################################################################################################


        ## Model flags
        self.Tbuildup_setting      = 1 # 1 = on # flag to control the buildup of turbulence along turbine rows, if disabled TI is equal to free stream TI
        self.wake_ind_setting      = 1 # 1 = on # flag to control the wake accumulation along turbine rows, if disabled U is equal to free stream velocity
        self.accu                  = 'dominant' # accumulation model
        self.accu_inlet            = True # accumulation of inlet velocity field. If enabled, the inflow flow field to each turbine is computed as the aggregated upstream flow field. If disabled, the model behaves like the HAWC2-DWM model [2] and [4]
        self.full_output           = False # control the amount of outputs scalar. If true, the velocity flow field is returned

        ## Misc variables
        self.mek_elec_power        = 0.95 # mek power to electric power
        self.rho                   = 1.25 # air at +10C
        self.wtg_ind               = np.asarray([]) # turbine index in main loop
        self.WTG_spec              = [] # class holding the turbine spec


        ## WF Control BEM
        self.derating              = 1.
        self.optim                 = False


        ## DWM Filter functions
        #if self.Keck:
        # Based on calibration data according to Keck et al. [5]. # F1 same as Madsen 2010
        self.f1                    =[0., 4.,] # F1 function = 1 after 4R
        self.f2                    =[0.035, 0.35] # F2 function starts at 0.035 for first 4R, then follows F2 = 1-0.965*e(-0.45*(2X-4)) (X in [R])
        """
        if self.Madsen:
            # Based on calibration data according to Keck et al. [5]. # F1 same as Madsen 2010
            self.f1 = [0., 4., ]  # F1 function = 1 after 4R

            # NOT the same F2!!!
            self.f2 = [0.035, 0.35]  # F2 function starts at 0.035 for first 4R, then follows F2 = 1-0.965*e(-0.45*(2X-4)) (X in [R])
        if self.Larsen:
            # Based on calibration data according to Keck et al. [5]. # F1 same as Madsen 2010
            self.f1 = [0., 4., ]  # F1 function = 1 after 4R

            # NOT the same F2!!!
            self.f2 = [0.035, 0.35]  # F2 function starts at 0.035 for first 4R, then follows F2 = 1-0.965*e(-0.45*(2X-4)) (X in [R])
        #"""

        # Iterative Progress for the Main Loop
        self.iT = 0.

    def parse(self,**par):
        """Parsing data to class holding all meta data for sDWM core model."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                print 'key: ', key
                raise Exception('Key is not found, check spelling')

        # --------------------------------------------------------------------------
        # ---------------------------------- # Model Specification Setting # ----------------------------------------- #

        if self.previous_sDWM:
            # Settings for Original sDWM
            self.steadyBEM_AINSLIE = False
            self.Meandering = True

            self.use_saved_data = False
            self.working_with_meandering_statistical_data = True
            self.previous_sDWM_working_with_a_MannBox = False  # we not use the original Meand matrix but data

            self.Keck = True

            ###################
            self.Madsen = False
            self.Larsen = False

            self.without_filter_functions = False
        if not self.previous_sDWM:
            self.steadyBEM_AINSLIE = True  # if True, BEM_Ainslie use the average deficit, False is not Implemented for now.

            # ----------- # Meandering # -------------------#
            self.Meandering = True
            # DO AN OPTION to DISABLE MEANDERING
            self.use_saved_data = True  # Use WAKES pre-computed and stored in the root 'src'
            self.working_with_meandering_statistical_data = False

            # ----------- # Wake Added Turbulence Settings # ----------- #
            self.WaT = True  # Traditional model of WaT based on  an isotropic Mannbox
            # In case of Keck model the WaT is managed in the deficit process, so it's not dynamic WaT


            # ------------- # Eddy Viscosity Model # ------------------ #
            self.EV_model = 'Larsen'
            # Keck : F1,F2 Improved WaT, Turbulent formulation, ABL etc..
            # Larsen: F1, F2, Famb, Ainslie
            # Madsen: F1, F2, Anslie or without F1, F2
            if self.EV_model == 'Keck':
                self.Keck = True
                self.Larsen = False
                self.Madsen = False

                self.without_filter_functions = False
            if self.EV_model == 'Larsen':
                self.Keck = False
                self.Larsen = True
                self.Madsen = False

                self.without_filter_functions = False
            if self.EV_model == 'Madsen':
                self.Keck = False  # F1,F2 Improved WaT, Turbulent formulation, ABL etc...
                self.Larsen = False  # F1, F2, Famb, Ainslie
                self.Madsen = True  # F1, F2, Anslie or without F1, F2

                self.without_filter_functions = False  # put True for 2008 Model
                # In the dynamic approach, we can average the deficits (and turbulence) in time,
                # and use this mean for the BEM-Ainslie Computation.
                # That's a drastic reduction of computation time.
                # However we can ask if it is really relevant for a dynamic simulation with a dynamic Input for Ainslie/BEM?


        ###### Model presented by Madsen in 2008 and 2010 Calibration for aeroelastic computation
        if self.Madsen:
            if self.without_filter_functions:
                # low turbulence intensity
                ###### Model presented by Madsen in 2008 Wake deficit and Turbulence simulated with two models
                # This model can be used to check kmt_r trend
                # That's a ref for the Wake added turbulence used in 2010 Madsen coalibrated model
                self.k_amb_Madsen = 0.001  # Madsen eddy viscosity without F1  (2008)
                self.k2_Madsen = 0.002  # Without F2
                self.km1 = 0.6
                self.km2 = 0.25
                self.kmt_r = []
            else:
                ###### Model presented by Madsen in 2010 Calibration for aeroelastic computation
                self.k_amb_Madsen = 0.07  # Madsen eddy viscosity with F1
                # self.k_amb_Madsen = 0.  # Madsen eddy viscosity with F1
                self.k2_Madsen = 0.008  # with F2

                def F2(X):
                    a1 = 0.046
                    b1 = -0.022

                    a2 = 0.821 / 25
                    # a2 = 1.
                    b2 = 0.0418 * 0.79
                    c2 = 0.2 - np.exp(0) + 0.6 + .3 + 0.1 / 10 * 3 + 0.005
                    Y = []
                    for x in X:
                        if x <= 2.:
                            y = 0.07
                        elif 2. < x < 7.:
                            y = a1 * x + b1
                        elif 7. <= x < 10.:
                            y = a2 * np.exp(b2 * x ** 2.) + c2
                            # y=0
                        elif x > 10.:
                            # y = 1.
                            y = 1.
                        Y.append(y)
                    return np.array(Y)

                self.F2 = F2
                # Same F1 that F1 used in Keck model
                # Not the same F2!

                self.km1 = 0.6
                self.km2 = 0.35
                self.kmt_r = []

        ###### Model presented by Larsen in 2012 Full scale validation of the DWM for loads and power production
        if self.Larsen:

            # Famb, F1 introduced/modified from 2010 to give better correlation with for more turbulence intensities cases
            # to model correctly the recovery
            self.k_amb_Larsen = 0.1
            self.k2_Madsen = 0.008  # seems to be the same as in 2010 (Calibrated Madsen Model)
            self.F_amb = 0.  # calculated in calc_mixl with 0.12/TI_amb

            self.km1 = 0.6
            self.km2 = 0.35
            self.kmt_r = []

            # raise Exception('Method not Implemented, no Formula for Famb')

            def F2(X):
                a1 = 0.046
                b1 = -0.022

                a2 = 0.821 / 25
                # a2 = 1.
                b2 = 0.0418 * 0.79
                c2 = 0.2 - np.exp(0) + 0.6 + .3 + 0.1 / 10 * 3 + 0.005
                Y = []
                for x in X:
                    if x <= 2.:
                        y = 0.07
                    elif 2. < x < 7.:
                        y = a1 * x + b1
                    elif 7. <= x < 10.:
                        y = a2 * np.exp(b2 * x ** 2.) + c2
                        # y=0
                    elif x > 10.:
                        # y = 1.
                        y = 1.
                    Y.append(y)
                return np.array(Y)

            self.F2 = F2

            # NEW F1 and Famb
            # I don't found the formula for f_amb


        # here we update the attributes that are affected by the new variable parsing
        self.mean_WS=[]
        self.mean_TI=[]
        self.mean_WS_DWM=[]
        self.mean_WS_rot=[]
        self.mean_TI_DWM=[]
        self.mean_TI_rot=[]
        self.power_estimate=[]
        self.bem_power_estimate=[]

        # self.hub_x                 = self.lx/2. # position of the rotor center
        self.hub_y                 = self.ly/2.
        self.nx                    = int(self.lx * self.dx)  # nb points in x (lateral)  direction global flow field
        self.ny                    = int(self.ly * self.dy)  # nb points in x (lateral)  direction global flow field
        self.nz                    = np.asarray([])  # nb points in z (streamwise) direction global flow field
        #print 'nx, ny, nz:'
        #print self.nx
        #print self.ny
        #print self.nz
        #raw_input(self.dx)
        self.x_vec                 = np.arange(1,self.nx+1,1)*(1.0/self.dx) # grid coordinates in lat. direction [R]
        self.y_vec                 = np.arange(1,self.ny+1,1)*(1.0/self.dy) # grid coordinates in long. direction [R]
        # Corrected to reach a real (1D*1D) domain in  case of perfectly aligned row
        #self.x_vec = np.linspace(0., self.lx,self.nx)# grid coordinates in lat. direction [R]
        #self.y_vec = np.linspace(0., self.ly, self.ny)  # grid coordinates in long. direction [R]
        #raw_input(self.x_vec)
        self.z_vec                 = np.asarray([]) # grid coordinates in streamwise. direction [D]
        self.z_vec_old             = np.asarray([]) # grid coordinates in streamwise. direction [D] at previous turbine
        self.x_mat                 = np.tile(self.x_vec.reshape(len(self.x_vec),1),self.ny)  # matrix grid holding x coordinates
        self.y_mat                 = np.tile(self.y_vec,(self.nx,1)) # matrix grid holding y coordinates
        self.z_mat                 = np.asarray([]) # matrix grid holding z coordinates

        #mixL Domain
        #we can increase self.lr_mixl
        self.lr_mixl               = self.lx+1. # 10 for single wake model, 6 for field model ? domain size in radial direction [R] (Should be about 0.5R wider flow field domain width)
        #print 'self.lr_mixl: ', self.lr_mixl
        self.vr_mixl               = np.linspace(0,self.lr_mixl-self.dr,(self.lr_mixl)/self.R_WTG*self.dr_mixl)  # coordinate vector
        self.nr_mixl = 0
        self.vr_m                  = np.arange(1,self.dr_mixl*self.lr_mixl+1,1,dtype=float) / self.dr_mixl  # coordinate vector for polar discretization
        self.dA_DWM                = np.concatenate(([0], (pi*self.vr_m[1:]**2 - pi*self.vr_m[0:-1]**2)), axis=0)  # area vector for polar discretization
        # self.nr                    = round(np.sqrt((max(self.x_vec - self.hub_x))**2 + (max(self.y_vec - self.hub_y))**2 ) * self.dR)
        # self.r_dist_2              = np.sqrt((self.x_mat - self.hub_x)**2 + (self.y_mat - self.hub_y)**2 )  #Find distance to centre of wake plane


class Meand:
     """ Class holding meandering magnitude related data
     """
     def __init__(self):
        """ Initializing the data / Set default value ."""

        ## Default Ambient conditions
        self.time                  = np.arange(1,600+1,1) # time vector in second for the meandering process
        self.x_offset              = np.asarray([]) # offset to wind dir relative row
        self.std_meand_x           = np.asarray([]) # std deviation lateral of wake center
        self.std_meand_y           = np.asarray([]) # std deviation long of wake center
        self.meand_pos_x           = np.asarray([]) # wake position vector lat
        self.meand_pos_y           = np.asarray([]) # wake position vector long
        self.rea                   = np.asarray([]) # fixed realization for debugging and unitest only

        self.WakesCentersLocations_in_time = []

        self.nt = len(self.time)

     def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')

class FFoR:
    # FFoR BOX is scaled per [D] /!\
     def __init__(self):
        """ Initializing the data / Set default value ."""
        self.x_vec_t               = np.asarray([])    #coordinate vector in fixed frame of reference nx,nz
        self.x_mat_t               = np.asarray([])    # coordinate matrix in fixed frame of reference (nx,ny,nz))

     def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')



class MFoR:
    """ Initialize structure holding calc results """

    def __init__(self,name):
        """ Initializing the data ."""
        self.name=name
        self.U=np.array([]) # streamwise velocity matrix in meandering frame of reference (vz_mixl,vr_mixl)
        self.V=np.array([] )# lateral velocity matrix in meandering frame of reference (vz_mixl,vr_mixl)
        self.visc=[]        # total viscosity matrix (vz_mixl,vr_mixl)
        self.TI_DWM=[]      # turbulence intensity formulated based on the relation derived between u'v' and u'u' (see Keck et al. [5])
        self.Turb_Stress_DWM=[] # matrix holding turbulence stress
        self.du_dr_tot=[]    # total velocity gradient (meandering + ABL)
        self.du_dr_DWM=[]    # velocity gradient from meandering
        self.Shear_add_du_dr=[] # added shear velocity gradient
        self.U_init =[]   # initial radial velocity deficit (axisymetric) input to ainslie model
        self.U_init_raw =[]
        self.U_init_pow =[]


    def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')


class Outputs:
    """ Initialize structure holding calc results """

    def __init__(self):
        """ Initializing the data ."""
        self.WS_axial_sym           = np.asarray([])
        self.WS_axial_sym_pow3      = np.asarray([])
        self.TI_axial_sym           = np.asarray([])
        self.TI_meand_axial_sym     = np.asarray([])
        self.TI_tot_axial_sym       = np.asarray([])
        # Store Mean turbine data
        self.C2C                     = np.asarray([])
        self.WS                     = np.asarray([])
        self.TI                     = np.asarray([])
        self.WS_DWM_BC              = np.asarray([])
        self.TI_DWM_BC              = np.asarray([])
        self.CP= np.asarray([])
        self.CPloc= np.asarray([])
        self.CQ= np.asarray([])
        self.CQlocCP= np.asarray([])
        self.CPloc= np.asarray([])
        self.CQ= np.asarray([])
        self.CQloc= np.asarray([])
        self.CT= np.asarray([])
        self.CTloc= np.asarray([])
        self.Cd= np.asarray([])
        self.Cl= np.asarray([])
        self.Cn= np.asarray([])
        self.Ct= np.asarray([])
        self.Edge= np.asarray([])
        self.F= np.asarray([])
        self.Flap= np.asarray([])
        self.Fperf= np.asarray([])
        self.Fshen= np.asarray([])
        self.Gamma= np.asarray([])
        self.PITCH= np.asarray([])
        self.PITCH_opt= np.asarray([])
        self.Pn= np.asarray([])
        self.Power= np.asarray([])
        self.Pt= np.asarray([])
        self.R= np.asarray([])
        self.RPM= np.asarray([])
        self.RPM_opt= np.asarray([])
        self.Re= np.asarray([])
        self.Set= np.asarray([])
        self.ThrLoc= np.asarray([])
        self.ThrLocLn= np.asarray([])
        self.Thrust= np.asarray([])
        self.Torque= np.asarray([])
        self.TqLoc= np.asarray([])
        self.TqLocLn= np.asarray([])
        self.Un= np.asarray([])
        self.Ut= np.asarray([])
        self.Vrel= np.asarray([])
        self.a= np.asarray([])
        self.a_last= np.asarray([])
        self.alpha= np.asarray([])
        self.aprime= np.asarray([])
        self.aprime_last= np.asarray([])
        self.nIt= np.asarray([])
        self.phi= np.asarray([])
        self.r= np.asarray([])
        self.uia= np.asarray([])
        self.uit= np.asarray([])

        self.U_init_raw=np.array([])
        self.U_init=np.array([])
        self.vr_mixl=np.array([])
        self.vz_mixl=np.array([])

        self.dA=np.array([])
        self.mean_a=np.array([])
        self.f_w=np.array([])
        self.r_w=np.array([])
        self.U_w=np.array([])

        self.Shear_add_du_dr=np.array([])
        self.Shear_add_du_dz=np.array([])
        self.TI_DWM=np.array([])
        self.Turb_Stress_DWM=np.array([])
        self.U=np.array([])
        self.U_init=np.array([])
        self.U_init_pow=np.array([])
        self.U_init_raw=np.array([])
        self.V=np.array([])
        self.WakeW=np.array([])
        self.du_dr_DWM=np.array([])
        self.du_dr_tot=np.array([])
        self.visc=np.array([])


    def parse(self,**par):
        """Parsing data ."""
        for key, value in par.iteritems():
            if hasattr(self,key):
                setattr(self, key, value)
            else:
                raise Exception('Key is not found, check spelling')


class Aero:

    """ Load inductions from available libraries based on turbine codename """
    def __init__(self,name):
        """ Initializing the data ."""
        self.name=name
        self.a=np.array([])
        self.r_t=np.array([])     # BEM radial discretization at rotor plane
        self.CP=np.array([])
        self.CPloc=np.array([]) # local CP on BEM grid
        self.Power=np.array([])
        self.pow_cur=np.array([])
        self.ct_cur=np.array([])
        self.dA=[]
        self.mean_a=[]
        self.f_w=[]
        self.r_w=[]
        self.U_w=[]
        self.Un=[]
        self.Ut=[]
        self.CT=np.array([])
        self.Cp=[]  # local Cp on flow field grid at rotor
        self.RPM=[]
        self.PITCH=[]
        self.RPM_opt=[]
        self.PITCH_opt=[]