# -*- coding: utf-8 -*-
""" Misc functions from variable sources
@moduleauthor:: Ewan Machefaux <ewan.machefaux@gmail.com>
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import pandas as pd

def LoadOutputs(folder,vWD,WF,WS,TI):
    powers=np.zeros((len(range(WF.nWT)),len(vWD)))
    ref_powers=[]
# for iD in range(0,1):
    for iD in np.arange(0,len(vWD),1):
        WD=float(vWD[iD])
        filename= folder +'/dwm_WS_'+ str(WS) +'_dir_' + str(WD)+'_TI_' + str(TI) +'.npy'
        # print filename
        try:
            out=np.load(filename).item()
            for iK in range(WF.nWT):
                # powers[iK,iD]=out[str(iK)][0] # BEM
                powers[iK,iD]=out[str(iK)][4] # power curve
            ref_powers.append(max(powers[:,iD]))
        except:
            for iK in range(WF.nWT):
                powers[iK,iD]=np.nan
            ref_powers.append(np.nan)

    return powers, ref_powers

# Pierre Elouan Rethore functions
def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

# Pierre Elouan Rethore functions
def ism360(v, endp):
    """Make sure the direction is in [0.0, 360]"""
    if np.isnan(v):
        return v
    if v>=0.0:
        if v>endp and v<360.0:
            return endp
        elif v>=360.0:
            return v - 360.0
        else:
            return v
    else:
        return ism360(v + 360.0, endp)

# Pierre Elouan Rethore functions
def my_rolling_deg(df, x='wd', y='eff', dwd=2.5):
    inte = interp1d(df[x], df[y])
    inte2 = lambda x_: inte(ism360(x_, df[x].max()))
    def filter_func(d):
        return {y:quad(inte2, d-dwd, d+dwd)[0]/(2.*dwd),x:d}
    return pd.DataFrame(map(filter_func, df[x]))


def find_intersections(A, B):
    # min, max and all for arrays
    amin = lambda x1, x2: np.where(x1<x2, x1, x2)
    amax = lambda x1, x2: np.where(x1>x2, x1, x2)
    aall = lambda abools: np.dstack(abools).all(axis=2)
    slope = lambda line: (lambda d: d[:,1]/d[:,0])(np.diff(line, axis=0))

    x11, x21 = np.meshgrid(A[:-1, 0], B[:-1, 0])
    x12, x22 = np.meshgrid(A[1:, 0], B[1:, 0])
    y11, y21 = np.meshgrid(A[:-1, 1], B[:-1, 1])
    y12, y22 = np.meshgrid(A[1:, 1], B[1:, 1])

    m1, m2 = np.meshgrid(slope(A), slope(B))
    m1inv, m2inv = 1/m1, 1/m2

    yi = (m1*(x21-x11-m2inv*y21) + y11)/(1 - m1*m2inv)
    xi = (yi - y21)*m2inv + x21

    xconds = (amin(x11, x12) < xi, xi <= amax(x11, x12),
              amin(x21, x22) < xi, xi <= amax(x21, x22) )
    yconds = (amin(y11, y12) < yi, yi <= amax(y11, y12),
              amin(y21, y22) < yi, yi <= amax(y21, y22) )

    return xi[aall(xconds)], yi[aall(yconds)]

def to_bool(value):
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(value).lower() in ("yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def unique2d(a):
    x, y = a.T
    b = x + y*1.0j
    idx = np.unique(b,return_index=True)[1]
    return a[idx]