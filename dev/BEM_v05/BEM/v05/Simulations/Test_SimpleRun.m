%% Initialization
InitClear;

PATH.DATA_WT='/home/ewan/Work/PostDoc/WindfarmControl/dev/BEM_v05/WT-data/';

PATH.BEM='/home/ewan/Work/PostDoc/DWMstandalone/BEM/Oldv05/BEM/';
PATH.WTlib='/home/ewan/Work/PostDoc/DWMstandalone/BEM/Oldv05/WTlib/';
PATH.Wind='/home/ewan/Work/PostDoc/DWMstandalone/BEM/Oldv05/Wind/';

require('WTlib','v06');
require('BEM','v05');
require('Wind','v01')
%% Params
nGrid=30;
% % sWT='NTK500p'; Format='hawc'; 
% sWT='NREL5MW';  Format='hawc'; 
% sWT='SB2XLL';  Format='xblade'; 

sWT='NY2'; Format='hawc'; 
%% Initialization
[ WT ]   = fInitWT( sWT , Format,PATH.DATA_WT);
[ WT ]   = fSetRotorGrid(nGrid,WT);
WT.Spec.vSIMRef
% lambda=3;

% U0=9.0;
U0=5.91359178053;

% pitch=interp1(WT.Spec.vSIMRef(:,1),WT.Spec.vSIMRef(:,3),U0,'cubic','extrap');
rpm=interp1(WT.Spec.vSIMRef(9:10,1),WT.Spec.vSIMRef(9:10,2),U0,'cubic','extrap');
% RPM=lambda*U0/WT.Rotor.R*60/2/pi;
% RPM=10.9135;

%%
[ Sim ]  = fInitSim( WT , [ U0  rpm  -1.5 ]  );
% [ Sim ]  = fInitSim( WT  );
[ Wind ] = fInitWind( Sim );

% Setting algorithm parameters
[ Algo ]   = fInitAlgo();
Algo.BEM.bTipLoss=1;
% Algo.bRough=0;
Algo.bReInterp=0;
Algo.bSwirl=1;
Algo.BEM.CTCorrection='GlauertCT';
Algo.BEM.SwirlMethod='Hawc';

%% Simulation
[ BEM ] = fRunBEM(WT,Sim,Wind,Algo);

%% Plotting results
% figure
% plot(BEM.r,BEM.Gamma)

% figure
% plot(BEM.r,BEM.Cl)
% 
% figure
% plot(BEM.r,BEM.a)
% 
% figure
% plot(BEM.r,BEM.phi)
%% compare

BEM.Power

% > BEM.Power
% 
% ans =
% 
%      3.934081227004653e+05

% dispatchFigs(1)

B = interp1([5.75 6.00]',[361.65363 410.90654 ]',U0,'cubic','extrap');
B
