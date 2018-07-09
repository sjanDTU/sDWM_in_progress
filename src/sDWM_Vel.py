__author__ = 'ewan'
import matplotlib as mpl
if mpl.get_backend<>'agg':
    mpl.use('agg')
# mpl.use('WxAgg')
# mpl.interactive(True)
import matplotlib.pyplot as plt
import numpy as np
# print '\nBackend:\n' + mpl.get_backend() + '\n'
import WindFarm as wf
import WindTurbine as wt
#import matplotlib.pyplot as plt
from DWM_flowfield_farm import DWM_main_field_model
from math import pi, sqrt, isnan
from DWM_GClarsenPicks import get_Rw
from DWM_init_dict import init
import time
from scipy import io, interpolate
from DWM_misc import smooth, to_bool
import matplotlib.pylab as plt
import matplotlib._cntr as cntr
import multiprocessing



###########################################################################

def sDWM(derating,kwargs,xind):
    # ttt = time.time()
#     WD,WS,TI,WTcoord,WTG,HH,R,stab,accum,optim,
    WD = kwargs.get('WD')
    WS = kwargs.get('WS')
    TI = kwargs.get('TI')
    WTcoord = kwargs.get('WTcoord')
    WTG = kwargs.get('WTG')
    HH = kwargs.get('HH')
    R = kwargs.get('R')
    stab = kwargs.get('stab')
    accum = kwargs.get('accum')
    optim = to_bool(kwargs.get('optim'))

    # WT = wt.WindTurbine('Vestas v80 2MW offshore','V80_2MW_offshore.dat',70,40)
    # WF = wf.WindFarm('Horns Rev 1','HR_coordinates.dat',WT)
    WT = wt.WindTurbine('Windturbine','../WT-data/'+WTG+'/'+WTG+'_PC.dat',HH,R)
    WF = wf.WindFarm('Windfarm',WTcoord,WT)

    if optim is True:
        print 'Performing optimization'
        WT.CP = np.load('../WT-data/'+WTG+'/'+WTG+'_CP.npy')
        WT.CT = np.load('../WT-data/'+WTG+'/'+WTG+'_CT.npy')
        WT.lambdaf3=WT.CP[:,:,0]
        WT.PITCH3=WT.CP[:,:,1]
        WT.CP3=WT.CP[:,:,2]
        WT.CT3=WT.CT[:,:,2]
        WT.CP3[WT.CP3>(16./27.)]=0
        WT.CP3[np.isnan(WT.CP3)]=0
        WT.CT3[np.isnan(WT.CT3)]=0
        WT.CPc = cntr.Cntr(WT.PITCH3,WT.lambdaf3,WT.CP3)
        WT.CTc=cntr.Cntr(WT.PITCH3,WT.lambdaf3,WT.CT3)
    elif optim is False:
        print 'Not Performing optimization'
        derating=1.0*np.ones((WF.nWT))
        WT.CP = None
        WT.CT = None

    # Scaling wind farm to NREL's rotor size
    if 'Lill' in WTcoord:
        WF.vectWTtoWT=WF.vectWTtoWT*(WT.R/46.5) # 46.5 is the Lillgrund nominal radius of SWT turbine


    # Compute distance WT to WT in mean flow coordinate system
    distFlowCoord, nDownstream, id0= WF.turbineDistance(WD)

    # Init dictionnaries
    deficits, turb, inlets_ffor, inlets_ffor_deficits,inlets_ffor_turb,out, DWM, ID_waked, ID_wake_adj, Farm_p_out, WT_p_out, Vel_out,WT_pitch_out,WT_RPM_out=init(WF)

    # Extreme wake to define WT's in each wake, including partial wakes
    ID_wake = {id0[i]:(get_Rw(x=distFlowCoord[0,id0[i],:],\
                              R=1.*WF.WT.R,TI=TI,CT=WT.get_CT(WS),pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0])>\
               np.abs(distFlowCoord[1,id0[i],:])).nonzero()[0] \
               for i in range(WF.nWT)}

    # Power output list
    Farm_p_out=0.
    WT_p_out=np.zeros((WF.nWT))
    WT_pitch_out=np.zeros((WF.nWT,2))
    WT_RPM_out=np.zeros((WF.nWT,2))
    Vel_out=np.zeros((WF.nWT))
    ## Main DWM loop over turbines
    for iT in range(WF.nWT):
        # Define flow case geometry
        cWT = id0[iT]
        #Radial coordinates in cWT for wake affected WT's
        x=distFlowCoord[0,cWT,ID_wake[cWT]]
        C2C   = distFlowCoord[1,cWT,ID_wake[cWT]]

        index_orig=np.argsort(x)
        x=np.sort(x)
        row= ID_wake[id0[iT]][index_orig]
        C2C=C2C[index_orig]

        # Wrapping the DWM core model with I/O
        par={
         'WS':WS,
         'TI':TI,
         'atmo_stab':stab,
         'WTG':'NREL5MW',      #'WTG':'NREL5MW',
         'WTG_spec': WT,
         'wtg_ind': row, # turbine index
         'hub_z':x/(2*WF.WT.R), # compute the flow field at downstream location of each downwind turbine !
         'hub_x': np.ceil((2*(max(abs(C2C))+WF.WT.R))/WF.WT.R)*0.5+C2C/(WF.WT.R), # lateral offset of each downwind turbine with respect to the most upstream turbine in the tree
         'C2C': C2C, # center 2 center distances between hubs
         'lx':np.ceil((2.*(max(abs(C2C))+WF.WT.R))/WF.WT.R), # length of the domain in lateral
         'ly':np.ceil((2.*(max(abs(C2C))+WF.WT.R))/WF.WT.R),  # length of the domain in longitudinal in D
         # 'wake_ind_setting':1,
         'accu_inlet': True,
         'derating': derating[row[0]],
         'optim': optim,
         'accu': accum, # type of wake accumulation
         'full_output': True # Flag for generating the complete output with velocity deficit
        }
        ID_wake_adj[str(id0[iT])]=row
        aero, meta, mfor, ffor, DWM, deficits,inlets_ffor,inlets_ffor_deficits, inlets_ffor_turb,turb, out,ID_waked = DWM_main_field_model(ID_waked,deficits,inlets_ffor,inlets_ffor_deficits,inlets_ffor_turb,turb,DWM,out,**par)
        # Farm_p_out= Farm_p_out+out[str(meta.wtg_ind[0])][4] # based on power curve

        # Total power
        Farm_p_out= Farm_p_out+out[str(meta.wtg_ind[0])][0] # based on BEM
        # Power by each turbine
        WT_p_out[iT]=out[str(meta.wtg_ind[0])][0]
        # Pitch and RPM
        WT_pitch_out[iT,0]=aero.PITCH
        WT_pitch_out[iT,1]=aero.PITCH_opt
        WT_RPM_out[iT,0]=aero.RPM
        print 'TSR = ' +str(WT_RPM_out[iT,0]*pi/30*63/15)
        WT_RPM_out[iT,1]=aero.RPM_opt
        Vel_out[iT]=out[str(meta.wtg_ind[0])][1]
        Velocity=ffor.WS_axial_ffor
        x=ffor.x_vec
        y=ffor.y_vec
    #
    #
    # print id0
    # print 'xind', xind
    # print 'wtp', WT_p_out[id0]
    # print 'pitch', WT_pitch_out[id0]
    # print 'omega', WT_RPM_out[id0]
    # print 'vel', Vel_out[id0]
    # print 'xind', xind[id0]


    print 'The farm production is: %4.2f kW, where each turbine is: %s' %(Farm_p_out,np.array_str(WT_p_out))
    return Velocity, x, y, meta