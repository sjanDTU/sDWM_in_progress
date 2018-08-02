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
    inputs = {'WD': 220,#270-70, #222
              'WS': 8.,
              'TI': 0.01,
              'WTcoord': '../data/Lill_rowB.dat',
              'WTG': 'NREL5MW',
              'HH': 65.0,
              'R': 46.5,
              'stab': 'N',
              'accum': 'dominant',
              'optim': False,
              'dynamic': True,
              'Meandering_turb_box_name': None,
              'WaT_turb_box_name': '1101'}
    ##########################################################################

    WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
    WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)


    P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])
