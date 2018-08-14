import numpy as np
import scipy.io as io
import csv #to save results in csv files

import WindFarm as wf
import WindTurbine as wt

from sDWM import sDWM
#from sDWM_stoch import sDWM_stoch




if __name__ == "__main__":
    inputs = dict()
    ##########################################################################
    inputs = {'WD': 0.,# ref is 270 in general (222 for lillgrund Row)  # Ref LES 180?or 0?
              'WS': 8.,
              'TI': 0.03,
              #'WTcoord': '../data/Lill_rowB.dat',
              #'WTcoord': '../data/TEST_S372m_nWT4_D0.dat',
              #'WTcoord': '../data/LES_WF_layout.dat',
              #'WTcoord': '../data/LES_WF_layout_3reduced.dat',
              'WTcoord': '../data/LES_WF_layout_7reduced.dat',
              #'WTcoord': '../data/LES_WF_layout_3reduced_modified.dat',
              'WTG': 'NREL5MW',
              'HH': 80.0,
              'R': 40.,
              'stab': 'N',
              'accum': 'dominant',
              'optim': False,
              'dynamic': True,                   # Imply bemsteady... TRUE
              #'Meandering_turb_box_name': ('1028','Mann_Box'),  # Imply used_saved data TRUE (not automatical for the moment)
              #'Meandering_turb_box_name': ('LES_Box_test','LES_Box'),
              'Meandering_turb_box_name': ('LES_Box_freq_5', 'LES_Box'),
              #'Meandering_turb_box_name': None,
              'WaT_turb_box_name': 'iso'}
    ##########################################################################

    WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
    WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)


    P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])
