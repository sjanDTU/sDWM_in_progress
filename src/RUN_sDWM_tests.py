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
              'WS': 9.,
              'TI': 0.06,
              #'WTcoord': '../data/Lill_rowB.dat',
              #'WTcoord': '../data/TEST_S372m_nWT4_D0.dat',
              #'WTcoord': '../data/LES_WF_layout.dat',
              'WTcoord': '../data/LES_WF_layout_3reduced.dat',
              'WTG': 'NY2',
              'HH': 80.0,
              'R': 40.,
              'stab': 'N',
              'accum': 'dominant',
              'optim': False,
              'dynamic': True,                   # Imply bemsteady... TRUE
              'Meandering_turb_box_name': None,  # Imply used_saved data TRUE (not automatical for the moment)
              'WaT_turb_box_name': '1028'}
    ##########################################################################

    WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
    WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)


    P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])
