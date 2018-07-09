import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io
import csv #to save results in csv files

import WindFarm as wf
import WindTurbine as wt

from sDWM import sDWM
#from sDWM_stoch import sDWM_stoch


#Stab Input
Stab=['VU','U','NU','N','NS']#,'S','VS']
#Wind Input
V_cut_in=6
V_cut_out=20

#Output for plot

P_TOT_STAB=[]
P_TOT=[]
P_IND=[]
for stab in Stab:
    P_TOT=[]
    for Wind_S in range(V_cut_in,V_cut_out+1):
        inputs = dict()
        ##########################################################################
        inputs = {'WD': 222,#270-70, #222
                  'WS': Wind_S,
                  'TI': 0.056,
                  'WTcoord': '../data/Lill_rowB.dat',
                  'WTG': 'NY2',
                  'HH': 65.0,
                  'R': 46.5,
                  'stab': stab,
                  'accum': 'dominant',
                  'optim': True}
        ##########################################################################

        WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
        WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)


        P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])

        # computation for plotting data

        P_TOT=P_TOT+[P_tot]
        P_IND=P_IND+[list(P_ind)]
    P_TOT_STAB=P_TOT_STAB+[P_TOT]


""" PLOT PART """
plt.figure(1)
for i in range(len(P_TOT_STAB)):
    #plt.subplot(3,3,i+1)
    #Plot Wind farm Power
    plt.plot(range(V_cut_in,V_cut_out+1), P_TOT_STAB[i], label=Stab[i]), plt.title(" Atmosphere Stability")
    plt.xlabel("Wind Speed (m/s)"),plt.ylabel("Power Production (kW)")
    #plt.show()
plt.legend()
plt.show()
"""
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
plt.plot(TI_Input,P_ind7)
plt.plot(TI_Input,P_ind6)
plt.plot(TI_Input,P_ind5)
plt.plot(TI_Input,P_ind4)
plt.plot(TI_Input,P_ind3)
plt.plot(TI_Input,P_ind2)
plt.plot(TI_Input,P_ind1)
plt.plot(TI_Input,P_ind0)
plt.title("Power production for each Turbine function of Turbulent Intensity")
plt.xlabel("TI")
plt.ylabel("Power (kW)")
plt.show()
"""







