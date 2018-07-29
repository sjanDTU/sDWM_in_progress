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

        # Discretization to read
        # Permit to compute the algo fastly with less points
        # to read with all points, discr_reduc_factor = 1, and nx_r = nx = 16384 points  2^14
        self.discr_reduc_factor = 1
        self.nx_r = 16384  # Discretization along x-axis to read fastly!

        self.SimulationTime = 30# (s)
        self.CorrectionDelay = True  # We want to begin the simulation when the first plan go out of the WindFarm Box
                                     # In this case at t=0, we have all the windfarm box affected by the turbulent box

        self.WakeExpansion = True

        self.Keck_Transport_Velocity = False  # True to apply 0.8*U advection transport: Keck synthetic turbulence


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
        self.RESULT_plot = True




    def Set(self, **parRotor):
        """Parsing data"""
        for key, value in parRotor.iteritems():
            setattr(self, key, value)