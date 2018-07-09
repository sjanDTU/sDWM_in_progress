import numpy as np
import scipy.io as io

import WindFarm as wf
import WindTurbine as wt

from sDWM import sDWM
from sDWM_stoch import sDWM_stoch



if __name__ == "__main__":
    inputs = dict()
    ##########################################################################
    inputs = {'WD': 222,
              'WS': 9,
              'TI': 0.056,
              'WTcoord': '../data/Lill_rowB.dat',
              'WTG': 'NY2',
              'HH': 65.0,
              'R': 46.5,
              'stab': 'N',
              'accum': 'dominant',
              'optim': True}
    ##########################################################################



    WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
    WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)

    dOpt = np.array([0.273399238 , 0.23551836 , 0.26183521 , 0.20363146 , 0.24715548 , 0.22095096 , 0.21827707 , 0.58795447])
    dOpt_relAvail = np.zeros(dOpt.shape)

    P_tot_opt, P_ind_opt, Pitch, RPM, Vel_out, id0 = sDWM(dOpt, inputs, [])

    dmax_it = np.ones((WF.nWT),)
    P_tot_norm, P_ind, Pitch, RPM, Vel_out, id0 = sDWM(dmax_it, inputs, [])
    dOpt_relAvail[id0[0]] = P_ind_opt[id0[0]]/ P_ind[id0[0]]

    for i in range(0,dOpt.shape[0]-1,1):
        dmax_it[id0[i]] =  dOpt[id0[i]]
        P_tot, P_ind,Pitch,RPM,Vel_out, id0 = sDWM(dmax_it,inputs,[])
        dOpt_relAvail[id0[i+1]] = P_ind_opt[id0[i+1]] / P_ind[id0[i+1]]

    print 'Optimum derating relative to turbine available power'
    print dOpt_relAvail

    print 'Power increase: ', P_tot_opt/ P_tot_norm

    # results = dict()
    # results = { 'Ptots': P_tot,
    #  'Pinds': P_ind,
    #  'Pitchs': Pitch,
    #  'RPMs': RPM,
    #  'Vels_out': Vel_out,}

    # io.savemat('../data/LillgrundEfficiency/sDWM/Rowwise/sDWM_rowB_WD222_WS9_TI56',results)