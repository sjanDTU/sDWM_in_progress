"""
Wed 23/05/2018, 15:03
Author:
Augustin Grosdidier
"""

import numpy as np

########################################################################################################################
def read_npy_WT_file_bem(WT_char,filename):
    Data = np.load('../WT-data/'+WT_char+'/BEM/'+filename+'.npy')
    return Data
########################################################################################################################
WT_char='NREL5MW'

#filename='pitch_opt'
#filename='PITCH'
#filename='Omega'
#filename='WScurve'
filename='Vsimref'

result = read_npy_WT_file_bem(WT_char,filename)
print result
print len(result)
print result[0]
print len(result[0])


########################################################################################################################
def read_npy_WT_file(WT_char,filename):
    Data = np.load('../WT-data/'+WT_char+'/'+filename+'.dat')
    return Data
########################################################################################################################
WT_char='NREL5MW'

#filename='pitch_opt'
#filename='PITCH'
#filename='Omega'
#filename='WScurve'
filename=WT_char+'_Spec'

result = read_npy_WT_file(WT_char,filename)
print result
print len(result)
print result[0]
print len(result[0])