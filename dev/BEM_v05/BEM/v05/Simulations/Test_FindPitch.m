%%
InitClear
require('BEM','v05');
require('WTlib','v06');
require('Wind','v01')
% sWT='SB2'; Format='xblade'; 
% sWT='Riso10MW'; Format='hawc'; 
% sWT='NREL5MW'; Format='hawc'; 
sWT='NY2'; Format='hawc'; 
% sWT='NTK500'; Format='hawc';
[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
WT.Nacelle.tilt=0;
%%
[ WT ]   = fSetRotorGrid(30,WT);
[ Algo ] = fInitAlgo();
Algo.bReInterp=0;

R=fWTFindPitch(WT ,4:0.5:25,WT.Spec.P_rated,[],1,0,0,[],Algo);
WT.Spec.vSIMRef=R.SIMRef;
fWTPowerCurve('BEM',WT ,4:0.1:25,1,Algo)

