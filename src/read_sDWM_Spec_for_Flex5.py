# -*- coding: utf-8 -*-
"""
Extracted from sDWM codes:

Created on Tue May  5 08:44:05 2015

@author: Ewan Machefaux
# Python interpreation of Matlab-based steady Blade Element Momentum of Emmanuel Branlard (emmanuel.branlard@gmail.com)
"""


import numpy as np
from math import pi





def fReadSpec(SpecFile):
    """Function that loads the turbine spec file
    """
    fd = open(SpecFile, 'r')
    WsRpmPitch=[]
    A=[]
    while 1:
        line = fd.readline()
        if not line:
            break
        pass
        if 'SPECS' in line:
            # Read RPM and Pitch curve
            while 1:
                line = fd.readline()
                if not line:
                    break
                else:
                    n = line.split(" ");
                    n = filter(None, n)
                    nn = [float(j) for j in n[0:3]]
                    WsRpmPitch.append(nn)
        else: # read other specs
            n = line.split(" ");
            n = filter(None, n)
            A.append(n[0])
    fd.close()
    WsRpmPitch=np.array(WsRpmPitch).T
    A=np.array(A)
    A = [float(j) for j in A]
    parRotor = {'nB': A[0],
                'BladeLength':A[1],
                'rhub': A[2],
                'cone':A[3],
                'R': A[1] + A[2]}

    parWT = {'tilt': A[4],
             'yaw': A[5],
             'H': A[6]}

    parAlgo = {'Ngrid': A[7]}

    parSpec = {'Omega_rated': A[8] * 2. * pi / 60.,
               'P_rated': A[9],
               'Pitch': WsRpmPitch[2,:],
               'RPM': WsRpmPitch[1,:],
               'WS': WsRpmPitch[0,:]}

    """
    PITCH = WsRpmPitch[2,:]
    RPM = WsRpmPitch[1,:]
    WS = WsRpmPitch[0,:]

    plt.figure(), plt.title('Wind Turbine Data'), plt.plot(WS,PITCH,label='Pitch(U)'), plt.legend(),plt.show()
    plt.figure(), plt.title('Wind Turbine Data'), plt.plot(WS, RPM, label='Omega(U)'), plt.legend(), plt.show()
    """
    return WsRpmPitch
#########################################################################################
#Spec = 'C:/Users/augus/OneDrive/Documents/Stage/Codes/sDWM/sDWM-WFC_Working_PyCharm_2.7.14/WT-data/NREL5MW/NREL5MW_Spec.dat'
#L = fReadSpec(Spec)
#print L
#########################################################################################
def convert_spec_to_data_for_Infile(WsRpmPitch):
    """
    (U(m/s),RPM(tr/mn),Pitch(rad))->(U(m/s), Pitch(deg),Pitch(rad))
    :param WsRpmPitch:
    :return:
    """
    new_WsPitch = list(WsRpmPitch)
    new_WsPitch[1] = WsRpmPitch[2]*180/pi
    return new_WsPitch

def add_Power_for_Infile_based_on_sDWM_pitch_law(WsPitch):
    """
    (U(m/s), Pitch(deg),Pitch(rad)) -> (U(m/s), Pitch(deg), Power(kW) from sDWM calculation )
    Run sDWM for 1st turbine to get Power
    :return:
    """
    WIND_S = WsPitch[0]

    result = Main_result_Wind(WIND_S)
    WsPitchPower = list(WsPitch)
    WsPitchPower[2] = result[:,1]

    print 'WIND_S', WIND_S
    print 'result[1]', result[1]

    return WsPitchPower

def convert_Array(L):
    lines = len(L[0])
    rows = len(L)
    print 'L', L
    M = np.zeros((lines, rows))
    print 'M',  M
    for i in range(rows):

        print 'M[:,i]', M[:, i]
        print 'L[i]', L[i]
        M[:,i] = L[i]
    return M

def Convert_Spec_to_Infile(WTG):
    Spec = 'C:/Users/augus/OneDrive/Documents/Stage/Codes/sDWM/sDWM-WFC_Working_PyCharm_2.7.14/WT-data/'+WTG+'/'+WTG+'_Spec_restrictedWind.dat'

    WsRpmPitch = fReadSpec(Spec)
    WsPitch = convert_spec_to_data_for_Infile(WsRpmPitch)
    WsPitchPower = add_Power_for_Infile_based_on_sDWM_pitch_law(WsPitch)
    WsPitchPower = convert_Array(WsPitchPower)
    np.savetxt(
        'C:/Users/augus/OneDrive/Documents/Stage/Codes/Comparisons_sDWMBEM_Flex5/WsPitchPower.dat',
        WsPitchPower)

    return WsPitchPower

##########################################################################################################################

# from Main_result_Wind import Main_result_Wind
#WTG = 'NREL5MW'
#WsPitchPower = Convert_Spec_to_Infile(WTG)
#print WsPitchPower