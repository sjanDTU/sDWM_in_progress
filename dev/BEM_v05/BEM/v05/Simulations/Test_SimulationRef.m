%%
clear all
close all
clc

tic
% InitClear
% setFigurePath('./')

global PATH
PATH.BEM='../../../../BEM_v05/BEM/';
PATH.WTlib='../../../../BEM_v05/WTlib/';
PATH.Wind='../../../../BEM_v05/Wind/';
% PATH.BEM='/home/ewan/Work/PostDoc/DWMstandalone/BEM/Oldv05/BEM/';
% PATH.WTlib='/home/ewan/Work/PostDoc/DWMstandalone/BEM/Oldv05/WTlib/';
% PATH.Wind='/home/ewan/Work/PostDoc/DWMstandalone/BEM/Oldv05/Wind/';
PATH.DATA_WT='../../../../BEM_v05/WT-data/';
require('BEM','v05');
require('WTlib','v06');
require('Wind','v01');

% sWT='SB2'; Format='xblade'; 
% sWT='Riso10MW'; Format='hawc'; 
% sWT='NREL5MW'; Format='hawc'; 
% sWT='NTK500'; Format='hawc'; 
% sWT='NY2'; Format='hawc'; 
sWT='V90'; Format='hawc'; 

[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
[ WT ]   = fSetRotorGrid(30,WT);
[ Algo ] = fInitAlgo();
Algo.BEM.bTipLoss=1;
Algo.bReInterp=0;
Algo.bSwirl=1;
Algo.BEM.CTCorrection='GlauertCT';
Algo.BEM.SwirlMethod='Hawc';
Algo.nbIt=200;
% WT.Rotor.cone=0;
WT.Nacelle.tilt=0;
WT.Rotor.cone=0;
%WT.Nacelle.tilt=0.0;
%WT.Rotor.cone=0.0;
%%
% --------------------------------------------------------------------------------
% --- Finding Turbine specifications 
% --------------------------------------------------------------------------------
Opts.TipSpeedMax=(WT.Rotor.R*pi*15.8)/30; % V90
% Opts.TipSpeedMax=(WT.Rotor.R*pi*15.947)/30; % Ny2
% Opts.TipSpeedMax=80; % NREL5MW
% Opts.OmegaMin=6.9*2*pi/60; % NREL5MW
% Opts.OmegaMin=7.6394*2*pi/60; % NY2
Opts.OmegaMin=12.5*2*pi/60; % V90
Opts.bOptimBelowPref=true;
% Opts.bOptimBelowPref=false;
Opts.WS_Startup=3.5;
% Opts.vWS_out=3:1:25;
% Opts.vLambda=1:1:20;
% Opts.vPitch=-12.:1:40.;

% Opts.vWS_out=3:1:25;
% Opts.vLambda=3:1:14;
% Opts.vPitch=-5.:1:25.;

Opts.vWS_out=3.5:1:25;
Opts.vLambda=2:1:9;
Opts.vPitch=-2.:1:27.;
Opts.Pref=WT.Spec.P_rated;
% Opts.Pref=1912300;


% Opts.TipSpeedMax=(WT.Rotor.R*pi*9.7009)/30; % NREL5MW
Opts.bPlot=1;
[R]= fWTSimRef(WT,Opts,Algo)

% Power=R.CPlambdaPitch.Power;
% RPM=R.CPlambdaPitch.RPM;
% Thrust=R.CPlambdaPitch.Thrust;
% lambdaf=R.CPlambdaPitch.lambda;
% PITCH=R.CPlambdaPitch.PITCH;
% Omega=R.CPlambdaPitch.Omega;
% CP=R.CPlambdaPitch.CP;
% CT=R.CPlambdaPitch.CT;

Power=R.CPlambdaPitch.Data.Power;
RPM=R.CPlambdaPitch.Data.RPM;
Thrust=R.CPlambdaPitch.Data.Thrust;
lambdaf=R.CPlambdaPitch.Data.lambda;
PITCH=R.CPlambdaPitch.Data.PITCH;
Omega=R.CPlambdaPitch.Data.Omega;
CP=R.CPlambdaPitch.Data.CP;
CT=R.CPlambdaPitch.Data.CT;


WScurve=R.WS;
Vsimref=R.vSIMRef;
CPmax=R.CPmax;
lambda_opt=R.lambda_opt;
OmegaMax=R.OmegaMax;
pitch_opt=R.pitch_opt;
Powcurve=R.Power;
RPMcurve=R.RPM;
CTcurve=[R.PowerCurveData.CT];

toc

WT.Spec.vSIMRef=R.vSIMRef;
R2=fWTPowerCurve('BEM',WT ,Opts.vWS_out,1,Algo)
% figure, fplotCodesComparison('WS','CT',{R2},'','','',1,1,[],[],'','');


%%
legds={}
%%
% disp('Watch out you are going to clear the workspace')
% pause
% load('WTSimRefWorkspace.mat')
% %%
% figure, fplotCodesComparison('r','CTloc',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),axis ij,box on,colorbar
% figure, fplotCodesComparison('r','a',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),axis ij,box on,colorbar
% figure, fplotCodesComparison('r','Gamma',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),box on,colorbar
%% Wind Farm control activities
% pitch shaving
% close all
% 
% figure()
% plot(R2.WS,R2.PITCH)
% 
% indo=find(R2.WS==13.5);
% R2.PITCH(indo);
% hold on
% pitch_Gradient=diff(R2.PITCH);
% indi=find(R2.WS==8.5);
% R2.PITCH(indi)
% % pitch_Gradient(ind)
% 
% plot(8.5,R2.PITCH(indi),'k+')
% plot(13.5,R2.PITCH(indo),'k+')
% 
% new_pitch=[R2.PITCH(1:indi) ;R2.PITCH(indo:end)];
% new_ws=[R2.WS(1:indi) ;R2.WS(indo:end)];
% 
% newpitchs=interp1(new_ws,new_pitch,R2.WS,'pchip');
% 
% plot(R2.WS,newpitchs,'g-')
% % new_pitch=smooth(R1PITCH,5)
% % plot(R2.WS,new_pitch)
% WT.Spec.vSIMRefold=WT.Spec.vSIMRef;
% figure()
% plot(R2.WS,R2.Thrust)

%% updating spec

% WT.Spec.vSIMRef(:,3)=newpitchs;
% 
% R3=fWTPowerCurve('BEM',WT ,Opts.vWS_out,1,Algo)
% figure, fplotCodesComparison('WS','CT',{R3},'','','',1,1,[],[],'','');


% test
% WT.Spec.vSIMRef=R.vSIMRef;
% vWS=[9.0]
% pitch=interp1(WT.Spec.vSIMRef(:,1), WT.Spec.vSIMRef(:,3),vWS,'cubic','extrap');
% rpm=interp1(WT.Spec.vSIMRef(:,1), WT.Spec.vSIMRef(:,2),vWS,'cubic','extrap');
% vSIM=[vWS(:) rpm(:) pitch(:)];
% [ Sim ]  = fInitSim(WT);
% [ Sim ]  = fSetSim( Sim, WT, vSIM ); 
% [ Wind ] = fInitWind(  ); 
% 
% 
% [ewan] = fBEMsteady(WT,Sim,Wind,Algo)

%% Calling Find pitch standalone
[ WT ]   = fInitWT( sWT, Format ,PATH.DATA_WT);
WT.Nacelle.tilt=5.;
WT.Rotor.cone=2.5;
[ WT ]   = fSetRotorGrid(30,WT);
[ Algo ] = fInitAlgo();
Algo.BEM.bTipLoss=1;
Algo.bReInterp=0;
Algo.bSwirl=1;
Algo.BEM.CTCorrection='GlauertCT';
Algo.BEM.SwirlMethod='Hawc';
Algo.nbIt=200;
Opts.TipSpeedMax=(WT.Rotor.R*pi*15.947)/30; % Ny2
% Opts.TipSpeedMax=80; % NREL5MW
% Opts.OmegaMin=6.9*2*pi/60; % NREL5MW
Opts.OmegaMin=7.6394*2*pi/60; % NY2
% Opts.bOptimBelowPref=true;
Opts.bOptimBelowPref=false;
% Opts.vWS=12;
% Opts.Pref=1912000;
% Opts.TipSpeedMax=(WT.Rotor.R*pi*9.7009)/30; % NREL5MW
% [R]= fWTFindPitch(WT,vSIMRef,Opts.vWS,Opts.Pref,0,Opts.bOptimBelowPref,Algo)
U=9; 
% Finding derated wind speed and RPM for given derated power
Pd   = 638100;
Pref = interp1(WScurve,Powcurve,U);
if Pd >= Pref
    warning('P derated above Prated')
end
Ud   = interp1(Powcurve,WScurve,Pd*0.001); % <<<<<<< POWER NOT UNIQUE
RPMd = interp1(WScurve,RPMcurve,Ud);
% Setting up the derated "WS-RPM" curve based on vSIMRef and RPM_derated%
vWSd=linspace(3,25,100);
vWSd=9;
[~,iUd]=min(abs(vWSd-Ud));
vRPM=interp1(WScurve,RPMcurve,vWSd); 
vRPM(iUd:end)=RPMd;

figure
plot(vWSd,vRPM)

vPITCHRef= zeros(1,length(vWSd)); % Unknown 
vSIMRef=[vWSd(:) vRPM(:) vPITCHRef(:)];

% R3=fWTFindPitch(WT, vSIMRef, Opts.vWS_out,Pd,0,Opts.bOptimBelowPref,Algo);
R3=fWTFindPitchEval(WT, vSIMRef, vWSd,Pd,0,Opts.bOptimBelowPref,Algo);
% [R]= fWTFindPitch(WT,vSIMRef,Opts.vWS,Opts.Pref,0,Opts.bOptimBelowPref,Algo)


%% Compt

figure();plot(R.WS,R.Power,'bo-')
hold on
plot(R2.WS,R2.Power,'r+-')
plot(R3.WS,R3.Power,'g+-')

figure();plot(R.WS,R.RPM,'bo')
hold on
plot(R2.WS,R2.RPM,'r+')
plot(R3.WS,R3.RPM,'g+-')

figure();plot(R.WS,R.PITCH,'bo-')
hold on
plot(R2.WS,R2.PITCH,'r+-')
plot(R3.WS,R3.PITCH,'g+-')

figure();plot(R.WS,[R.PowerCurveData.CT],'bo')
hold on
plot(R2.WS,R2.CT,'r+')
plot(R3.WS,R3.CT,'g+-')

figure()
plot(R2.CP,R2.CT,'r+')
max(R2.CP)
max(max(R.CPlambdaPitch.CP))
max(max(R.CPlambdaPitch.Data.CP))