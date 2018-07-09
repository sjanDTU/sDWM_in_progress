import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io
import csv #to save results in csv files

import WindFarm as wf
import WindTurbine as wt

from sDWM import sDWM
#from sDWM_stoch import sDWM_stoch


#TI Input
TI_cut_in=0.04
TI_cut_out=0.5
N=30
TI_Input=[]
for i in range(N):
    TI_Input=TI_Input+[TI_cut_in+i*(TI_cut_out-TI_cut_in)/float(N)]


#Output for plot

P_TOT=[]
P_IND=[]

for TI_input in TI_Input:
    inputs = dict()
    ##########################################################################
    inputs = {'WD': 222,#270-70, #222
              'WS': 10,
              'TI': TI_input,
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


    P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])

    # computation for plotting data

    P_TOT=P_TOT+[P_tot]
    P_IND=P_IND+[list(P_ind)]

""" PLOT PART """
plt.figure(1)
plt.subplot(1,2,1)
#Plot Wind farm Power
plt.plot(TI_Input, P_TOT), plt.title("Wind Farm Power Production function of Turbulent Intensity")
plt.xlabel("TI"),plt.ylabel("Power Production (kW)")
#plt.show()

# Power for each Turbines
P_ind0=[]
P_ind1=[]
P_ind2=[]
P_ind3=[]
P_ind4=[]
P_ind5=[]
P_ind6=[]
P_ind7=[]
for i in range(len(TI_Input)):

    P_ind7 = P_ind7 + [P_IND[i][0]]
    P_ind6 = P_ind6 + [P_IND[i][1]]
    P_ind5 = P_ind5 + [P_IND[i][2]]
    P_ind4 = P_ind4 + [P_IND[i][3]]
    P_ind3 = P_ind3 + [P_IND[i][4]]
    P_ind2 = P_ind2 + [P_IND[i][5]]
    P_ind1 = P_ind1 + [P_IND[i][6]]
    P_ind0 = P_ind0 + [P_IND[i][7]]

plt.subplot(1,2,2)
plt.plot(TI_Input,P_ind7,label='Turbine 0')
plt.plot(TI_Input,P_ind6,label='Turbine 1')
plt.plot(TI_Input,P_ind5,label='Turbine 2')
plt.plot(TI_Input,P_ind4,label='Turbine 3')
plt.plot(TI_Input,P_ind3,label='Turbine 4')
plt.plot(TI_Input,P_ind2,label='Turbine 5')
plt.plot(TI_Input,P_ind1,label='Turbine 6')
plt.plot(TI_Input,P_ind0,label='Turbine 7')
plt.title("Power production for each Turbine function of Turbulent Intensity")
plt.xlabel("TI"), plt.ylabel("Power (kW)")
plt.legend(loc=4)
plt.show()








