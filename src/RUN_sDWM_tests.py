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
    inputs = {'WD': 222,#270-70, #222
              'WS': 9.,
              'TI': 0.056,
              'WTcoord': '../data/Lill_rowB.dat',
              'WTG': 'NY2',
              'HH': 65.0,
              'R': 46.5,
              'stab': 'N',
              'accum': 'dominant',
              'optim': True,
              'dynamic': False}
    ##########################################################################

    WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
    WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)


    P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])
    #print Vel_out
    #print P_ind

    # P_tot_stoch, P_ind_stoch, Vel_out_stoch, xind = sDWM_stoch(np.ones((WF.nWT), ), inputs, [])

    #print 'Ptots', P_tot
    #print 'Ptots_stoch', P_tot_stoch
    #print 'Pinds', P_ind
    # print 'Pinds_stoch', P_ind_stoch
    #print 'Vels_out', Vel_out
    # print 'Vels_out_stoch', Vel_out_stoch
    #print 'Pitch', Pitch
    #print 'Pitch for first Turbine', Pitch[0,0]
    #print 'RPM', RPM

    """# Tests
    dmax = np.ones((WF.nWT))
    for iWT in np.arange(0., WF.nWT, 1.):
        dmax[iWT] = P_ind[iWT] / WT.P_rated

    print 'dmax:', dmax
    P_tot, P_ind, Pitch, RPM, Vel_out, id0 = sDWM(dmax, inputs, [])

    dOpt = np.array([ 0.26596013 , 0.23551836 , 0.26183521 , 0.20363146 , 0.24715548 , 0.22095096 , 0.21827707 , 0.58795447])

    P_tot_opt, P_ind, Pitch, RPM, Vel_out, id0 = sDWM(dOpt, inputs, [])

    print 'Ptots', P_tot
    print 'Pinds', P_ind

    print 'Power increase: ', P_tot_opt/ P_tot
"""
    """ 
    results = dict()
    results = { 'Ptots': P_tot,
    'Pinds': P_ind,
    'Pitchs': Pitch,
    'RPMs': RPM,
    'Vels_out': Vel_out,}

    io.savemat('../data/LillgrundEfficiency/sDWM/Rowwise/sDWM_rowB_WD222_WS9_TI56_test',results)
    #dict = {'Python': '.py', 'C++': '.cpp', 'Java': '.java'}
    fieldnames=[]
    for key in results:
        fieldnames=fieldnames+[key]
    import pandas as pd

    #pd.DataFrame(results).to_csv('OUTPUT.csv')

    #import csv

    with open('output.csv', 'wb') as output:
        writer = csv.writer(output)
        for key, value in results.iteritems():
            writer.writerow([key])
            #
            # value=[value]
            #print value
            if key=='Ptots':
                value=[value]
            for val in value:
                print val

                writer.writerow([val])"""