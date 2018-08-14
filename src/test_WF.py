import numpy as np
import matplotlib.pyplot as plt


"""
# Wake Added Turbulence Implementation

P_Madsen = [4963.7, 939.03, 1279.2,  1608.85, 1852.97, 1960.41, 1985.94, 2019.74]
P_Keck = [4963.7,  1422.91, 1528.75, 1665.02, 1721.63, 1684.75, 1713.17, 1736.01]
P_Larsen = [4963.7,  1115.26, 1072.04, 1242.05, 1491.85, 1685.4,  1822.07, 1958.84]

WTs = range(0,8)

plt.figure()
plt.plot(WTs, P_Keck, 'o', label='Keck')
plt.plot(WTs, P_Madsen, 'x', label = 'Madsen')
plt.plot(WTs, P_Larsen, '*', label = 'Larsen')
plt.xlabel('WTs'), plt.ylabel('P[kW]')
plt.legend()
plt.show()
#"""


L = np.loadtxt('../data/Lill_rowB.dat')
L = np.loadtxt('../data/row_test.dat')
L = np.loadtxt('../data/TEST_S186m_nWT16_D0.dat')
print L


"""
print np.arccos((L[:,0]-L[-1,0])/np.sqrt((L[:,0]-L[-1,0])**2+(L[:,1]-L[-1,1])**2))*180./3.1459

plt.figure()
plt.plot(L[:,0]-L[-1,0], L[:,1]-L[-1,1],'x')
plt.plot(L[:,0]-L[-1,0], L[:,1]-L[-1,1])
plt.axis('equal')
plt.show()
#"""
"""
L[:,0] = L[:,0]-L[-1,0]
L[:,1] = L[:,1]-L[-1,1]



np.savetxt('../data/row_test.dat', L)
#"""

#"""
import WindFarm as wf
import WindTurbine as wt
WTG = 'NY2'
HH = 90
R = 46.5

WD = 270

WTcoord = '../data/row_test.dat'
#WTcoord= '../data/Lill_rowB.dat'
#print 'RowB'

WT = wt.WindTurbine('Windturbine','../WT-data/'+WTG+'/'+WTG+'_PC.dat',HH,R)
WF = wf.WindFarm('Windfarm',WTcoord,WT)
print 'vecttotvecWT: ', WF.vectWTtoWT
distFlowCoord, nDownstream, id0= WF.turbineDistance(WD)

print 'distflowcoord', distFlowCoord
print 'ndownstream', nDownstream
print 'id0', id0
#"""

#Create a data file for WT coord. beginning with WT (0,0)
#According to a spacing, nb of WT, direction

def create_WT_row(spacing,nWT,D,row_name):
    """
    Create a data file for WT coord for a perfect row of WT. beginning with WT (0,0)
    According to a spacing, nb of WT, direction

    :param spacing: for the moment we admit a constant spacing (m)
    :param nWT: number of Turbine
    :param D: Direction of the row in degree (degree)
    :return:
    """
    L = [[0,0]]

    spacing_x = np.cos(D)*spacing
    spacing_y = np.sin(D)*spacing

    for i_WT in range(1,nWT):

        L.append([i_WT*spacing_x, i_WT*spacing_y])

    print L
    print row_name
    np.savetxt('../data/'+row_name+'.dat',L)
    return

#Test
if False:
    spacing = 46.5*4*2
    nWT = 6
    D = 0
    row_name = 'TEST_S'+str(int(spacing))+'m_nWT'+str(nWT)+'_D'+ str(int(D))

    create_WT_row(spacing, nWT, D, row_name)

