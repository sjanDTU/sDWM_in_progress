import numpy as np
import matplotlib.pyplot as plt
import scipy.io as io
import csv #to save results in csv files

import WindFarm as wf
import WindTurbine as wt

from sDWM import sDWM
#from sDWM_stoch import sDWM_stoch
#Wind Input


def Main_result_Wind(WIND_S):
    #Wind Input


    #Output for plot

    P_TOT=[]
    P_IND=[]
    PITCH_for_1st_turbine=[]
    RPM_for_1st_turbine=[]

    #"""



    #"""

    Wind_Power_Pitch = np.zeros((len(WIND_S),3))

    for Wind_S in WIND_S:
        inputs = dict()
        ##########################################################################
        inputs = {'WD': 222,#270-70, #222
                  'WS': Wind_S,
                  'TI': 0.1,
                  'WTcoord': '../data/Lill_rowB.dat',
                  'WTG': 'NREL5MW',#NY2 or NREL5MW
                  'HH': 90.,
                  'R': 63.,
                  'stab': 'N',
                  'accum': 'dominant',
                  'optim': False}
        ##########################################################################

        WT = wt.WindTurbine('Windturbine','../WT-data/'+inputs['WTG']+'/'+inputs['WTG']+'_PC.dat',inputs['HH'],inputs['R'])
        WF = wf.WindFarm('Windfarm',inputs['WTcoord'],WT)


        P_tot, P_ind,Pitch,RPM,Vel_out, id0= sDWM(np.ones((WF.nWT),),inputs,[])

        # computation for plotting data

        P_TOT=P_TOT+[P_tot]
        P_IND=P_IND+[list(P_ind)]
        PITCH_for_1st_turbine=PITCH_for_1st_turbine+[Pitch[7,0]]
        RPM_for_1st_turbine=RPM_for_1st_turbine+[RPM[0,0]]
        print RPM[0,0]


    """ PLOT PART """
    plt.figure(1)
    plt.subplot(1, 2, 1)



    #Plot Wind farm Power
    plt.plot(WIND_S, P_TOT), plt.title("Wind Farm Power Production function of Wind speed")
    plt.xlabel("Wind Speed (m/s)"),plt.ylabel("Power Production (kW)")
    #plt.show()

    #print P_IND
    #"""
    # Power for each Turbines
    P_ind0=[]
    P_ind1=[]
    P_ind2=[]
    P_ind3=[]
    P_ind4=[]
    P_ind5=[]
    P_ind6=[]
    P_ind7=[]
    for i in range(len(WIND_S)):

        P_ind7 = P_ind7 + [P_IND[i][0]]
        P_ind6 = P_ind6 + [P_IND[i][1]]
        P_ind5 = P_ind5 + [P_IND[i][2]]
        P_ind4 = P_ind4 + [P_IND[i][3]]
        P_ind3 = P_ind3 + [P_IND[i][4]]
        P_ind2 = P_ind2 + [P_IND[i][5]]
        P_ind1 = P_ind1 + [P_IND[i][6]]
        P_ind0 = P_ind0 + [P_IND[i][7]]

        Wind_Power_Pitch[i, 0] = WIND_S[i]
        Wind_Power_Pitch[i, 1] = P_ind0[i]
        Wind_Power_Pitch[i, 2] = PITCH_for_1st_turbine[i]
    plt.subplot(1, 2, 2)
    plt.plot(WIND_S,P_ind7,label='Turbine 0')
    plt.plot(WIND_S,P_ind6,label='Turbine 1')
    plt.plot(WIND_S,P_ind5,label='Turbine 2')
    plt.plot(WIND_S,P_ind4,label='Turbine 3')
    plt.plot(WIND_S,P_ind3,label='Turbine 4')
    plt.plot(WIND_S,P_ind2,label='Turbine 5')
    plt.plot(WIND_S,P_ind1,label='Turbine 6')
    plt.plot(WIND_S,P_ind0,label='Turbine 7')
    plt.title("Power production for each Turbine function of Wind Speed")
    plt.xlabel("Wind velocity (m/s)")
    plt.ylabel("Power (kW)")
    plt.legend(loc=4)
    plt.show()
    #"""
    plt.figure(),plt.plot(WIND_S,PITCH_for_1st_turbine,label='Pitch(U)'),plt.legend(),plt.show()

    np.savetxt(
        'C:/Users/augus/OneDrive/Documents/Stage/Codes/Comparison_BEM_sDWM/sDWM_Power_and_Pitch_for_first_Turbine_Flex5Pitch.dat',
        Wind_Power_Pitch)
    #np.savetxt('C:/Users/augus/OneDrive/Documents/Stage/Codes/Comparison_BEM_sDWM/sDWM_Power_and_Pitch_for_first_Turbine.dat', Wind_Power_Pitch)
    return Wind_Power_Pitch

"""
V_cut_in=8.
V_cut_out=16.
WIND_S=np.linspace(V_cut_in,V_cut_out,30)
#"""
#"""
WTG = 'NREL5MW'
INPUT = 'sweep'
dir = 'C:/Users/augus/OneDrive/Documents/Stage/Results/Flex5/'+WTG+'/'+INPUT+'/Data/'
kind = 'tim'

DATA = np.loadtxt(dir+'converged_'+kind+'.dat')  # (s), (m/s), (deg), (kW), (...)...
WIND_S = DATA[7:,1]
#"""
"""
from read_sDWM_Spec_for_Flex5 import fReadSpec
WTG='NREL5MW'
#########TOTAL Spec file########
#Spec = 'C:/Users/augus/OneDrive/Documents/Stage/Codes/sDWM/sDWM-WFC_Working_PyCharm_2.7.14/WT-data/'+WTG+'/'+WTG+'_Spec.dat'
######Restricted Spec file#########
Spec = 'C:/Users/augus/OneDrive/Documents/Stage/Codes/sDWM/sDWM-WFC_Working_PyCharm_2.7.14/WT-data/'+WTG+'/'+WTG+'_Spec_restrictedWind.dat'
L = fReadSpec(Spec)
WIND_S = L[0]
#"""
#"""
result = Main_result_Wind(WIND_S)
#"""
"""
plt.figure(1)
plt.subplot(121), plt.plot(result[0], result[1], label = 'Power(U)'), plt.legend()
plt.subplot(122), plt.plot(result[0], result[2], label = 'Pitch(U)'),plt.plot(L[0], L[2], label='Pitch(U) from Spec file'), plt.legend()
plt.show()
#"""