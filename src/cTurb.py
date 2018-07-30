"""
author: Augustin Grosdidier
date: 07/06/2018

Purpose:
Define some class variable to store data along the code computation
"""

class MannBox:

    def __init__(self):
        # Required Params
        self.fieldDim = 0.
        self.N_Comp = 0.

        self.nx = 1      # Point Discritization along x-axis
        self.ny = 1      # Point Discritization along y-axis
        self.nz = 1      # Point Discritization along z-axis


        self.lx = 1.     # 32768m in total
        self.ly = 1.
        self.lz = 1.

        self.L = 1.
        self.U = 1.
        self.T = 1.

        self.dx = 1.
        self.dt = 1.

        self.ti = []

        self.TurbBox_of_interest = []
        self.u_TurbBox = []  # For wake added Turbulence
        self.v_TurbBox = []
        self.w_TurbBox = []
        self.plan_of_interest = []
        self.radius_plan_of_interest = []
        self.values_of_interest= []
        #self.CenterLocation_on_Plan = 0., 0.
        self.WakeRadius_for_this_Plan = 0.
        self.Af = 0.


        self.U_ref = 0.
        self.R_ref = 0.
        self.TI = 0.

        self.TI_u = None

        # Plot option setting
        self.RESULT_plot = False

        # ------------------------------- # Method Settings # -------------------------------------------------------- #

        # Discretization to read
        # Permit to compute the algo fastly with less points
        # to read with all points, discr_reduc_factor = 1, and nx_r = nx = 16384 points  2^14
        self.discr_reduc_factor = 1
        self.nx_r = 16384  # Discretization along x-axis to read fastly!

        self.SimulationTime = 60  # (s)
        self.CorrectionDelay = True  # We want to begin the simulation when the first plan go out of the WindFarm Box
        # In this case at t=0, we have all the windfarm box affected by the turbulent box
        self.delay = 0.

        self.WakeExpansion = True  # Carry by a simple wake model proposed by Larsen 2009

        self.Keck_Transport_Velocity = False  # True to apply 0.8*U advection transport: Keck synthetic turbulence

        self.loop_debug = False
        self.multiplewake_build_on_first_wake = True

        self.Box_Kind = 'MannBox'  # LESBox or MannBox




    def Set(self, **parRotor):
        """Parsing data"""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)

class interpo_integrate():

    def __init__(self):
        # Required Params
        self.Interp_method = 'RBS'  # RBS is the fastest, but we have to implement a transpose for the good result
                                         # , or interp2d
        if self.Interp_method == 'interp2d':
            self.Interp2dKind = 'linear'
        self.f_cart = 0.                 # Interpolated function, assigned along the algo

        # Discretization for trapz on circle area, for the moment nr have to be equaled to nt
                                # 10-5   #10-3
        self.nr = 50            # 110      15
        self.nt = 50           # 100      15
        self.r = []
        self.t = []

        self.PolMesh = 0.
        self.CartMesh = 0.

        self.PolMesh_WakeCentered = 0.
        self.CartMesh_WakeCentered = 0.

        self.ValMesh = 0.

        # Interpolated function to estimate wakecenter at a specific distance, written in term of time
        self.f_vc_t = 0.
        self.f_wc_t = 0.
        self.F_tm_fvc_fwc = []

        #Plot Setting Options
        self.Plot = False


    def Set(self, **parRotor):
        """Parsing data ."""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)

class meand_mann():

    def __init__(self):
        print 'meand_mann'
        self.init_vc_wc = []
        self.yc_zc = []
        self.wake_center_location = (0, 0)
        self.ti_wake_center_locations = {}
        self.wakecenterlocation_in_time_ground_referential = []

    def Set(self, **parRotor):
        """Parsing data ."""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)

class windfarm():
    def __init__(self):
        self.U_mean = 0.

        self.WT_R = 0.      # List of radius, but i think for the moment we should work with just one size, so a float!

        self.WT_Rw = []     # for each WT, handle the list of Rw according to the distance in the flowcoord
                            # but for init we should work with a float!
        self.lenght = 0.

        self.TI = 0.
        self.CT = 0.   #One value for all the wind Farm

        self.stream_location_z = []
        self.nodim_lenght= 0.
        self.Rw_Max = 0.

        self.constant_spacing = True


    def Set(self, **parRotor):
        """Parsing data ."""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)
