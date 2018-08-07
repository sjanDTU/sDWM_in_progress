"""
author: Augustin Grosdidier
date: 02/08/2018

Purpose:
Create Data file for WT coord in a WindFarm.
This data file can be used instead of Lillgrund Data for a simulated WF

For the moment we just create a perfect row of Turbine (perfectly aligned), with same spacing.
"""
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
# Generate a WindFarm according to spacing, nWT, row direction
#"""
if True:
    spacing = 46.5*4*2
    nWT = 4
    D = 0
    row_name = 'TEST_S'+str(int(spacing))+'m_nWT'+str(nWT)+'_D'+ str(int(D))

    create_WT_row(spacing, nWT, D, row_name)

L = np.loadtxt('../data/Lill_rowB.dat')
L = np.loadtxt('../data/row_test.dat')
L = np.loadtxt('../data/TEST_S186m_nWT16_D0.dat')
print L
#"""

# To visualize the Wind Farm
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

# What is done in sDWM with the file
"""
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

# -------------------------- # Read Rot file to extract the layout of a windFarm # ----------------------------------- #

# In this case the layout of the LES windFarm
def plot_WF(XYZ_coord):
    print np.shape(XYZ_coord)
    plt.figure()
    plt.title('WindFarm Layout')
    plt.plot(XYZ_coord[:,0], XYZ_coord[:,2], '*')
    plt.xlabel('X position [R]'), plt.ylabel('Z position [R]')
    plt.show()
    return
def extract_coord_from_ROT_file(rot_file_name):
    path = '../LES_data/'
    f = open(path+rot_file_name, 'r')
    #print f.readlines()
    bool_r = False
    L = []
    for line in f:
        if '# Wake turbulence box' in line:
            bool_r = False

        if bool_r:
            L.append(line.split(' '))
            #print line

        if '# Row  Col     X        Y        Z' in line:
            bool_r = True

    #print L

    for i_l in range(len(L)):
        l_new = []
        for ll in L[i_l]:
            if ll !='' and ll!='\n':
                l_new.append(ll)
        #print l_new
        L[i_l] = l_new
    L = np.array(L)
    L = L.astype(float)

    XYZ_coord = L[:, 2:5]
    print XYZ_coord

    return XYZ_coord


# To generate a WindFarm according to a ROT file
#"""
ROT_name = 'v8z12_ABL_16wt_TURB3new_19.ROT'

XYZ_coord = extract_coord_from_ROT_file(ROT_name)
plot_WF(XYZ_coord)

R = 40. # KNOW THE R HERE

XZ_coord = XYZ_coord[:,(0,2)] * R
XZ_coord[:, 0] = 0.
np.savetxt('../data/LES_WF_layout.dat',XZ_coord)
print 'Coord Saved with name LES_WF_layout'
# X Lateral position
# Y Height Position
# Z Downstream Position

# Let's consider same Height for all WT
# To give the sDWM WTcoord we just need X, Z position

#"""