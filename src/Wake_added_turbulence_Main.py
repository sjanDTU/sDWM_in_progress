"""
date: 24/07/18
author: Augustin Grosdidier
Purpose: Treatment of a MannBox for wake added turbulence process
"""

import numpy as np
import math
import matplotlib.pyplot as plt

from cMann import MannBox

from ReadMannTurbulence import *

def pre_init_turb_WaT(filename):
    """
    Purpose:
    Extract Ti_u for wake added Components
    :return:
    """

    # ------------------ # GET TURBULENT BOX # -------------------------- #
    MannBox = ReadMannInput(filename)

    u_TurbBox = get_turb_component_from_MannBox(filename, 'ufluct', False, MannBox, video=False)
    U_mean = get_averaged_U(filename)# Specified U_mean for this Mannbox
    TI_u = np.sqrt(np.mean((u_TurbBox/U_mean) ** 2))


    return MannBox

# Test
result = pre_init_turb_WaT('1028')
print result.TI_u
