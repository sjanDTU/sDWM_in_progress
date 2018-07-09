function [ Power, RPM, Thrust, lambdaf, PITCH, Omega, CP, CT, WScurve, Vsimref, CPmax, lambda_opt, OmegaMax, pitch_opt, Powcurve, RPMcurve,CTcurve, CPcurve,Alpha ] = fSimulationRef()
%FSIMULATIONREF Summary of this function goes here
%   Detailed explanation goes here

% %%
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
sWT='NREL5MW'; Format='hawc'; 
% sWT='NTK500'; Format='hawc'; 
% sWT='NY2'; Format='hawc'; 

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
WT.Nacelle.tilt=5.;  % NREL 5MW
WT.Rotor.cone=2.5; % NREL 5MW
% WT.Nacelle.tilt=5.;  % NY2
% WT.Rotor.cone=2.7; % NY2
%WT.Nacelle.tilt=0.0;
%WT.Rotor.cone=0.0;
%%
% --------------------------------------------------------------------------------
% --- Finding Turbine specifications 
% --------------------------------------------------------------------------------
% Opts.TipSpeedMax=(WT.Rotor.R*pi*15.947)/30; % Ny2
Opts.TipSpeedMax=80; % NREL5MW
Opts.OmegaMin=6.9*2*pi/60; % NREL5MW
% Opts.OmegaMin=7.6394*2*pi/60; % NY2
Opts.bOptimBelowPref=true;
% Opts.bOptimBelowPref=false;% NREL5MW
Opts.WS_Startup=3.0;% NREL5MW
Opts.vWS_out=3:0.5:25;% NREL5MW
% Opts.WS_Startup=4.0;
% Opts.vWS_out=4:0.25:25;
Opts.vLambda=0.5:0.5:20;
Opts.vPitch=-5.:0.5:35.;
% Opts.vWS_out=3:2:25;15
% Opts.vLambda=2:2:15;
% Opts.vPitch=-15.:2:35.;
Opts.Pref=WT.Spec.P_rated;
Opts.bPlot=0;
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
CPcurve=[R.PowerCurveData.CP];


M=[R.CPlambdaPitch.Data.Results];
for i=1:length(Opts.vPitch); for j=1:length(Opts.vLambda); Alpha(i,j)=mean(M(i,j).alpha(30/2:30));end;end


% --------------------------------------------------------------------------------
% --- Computing Power curve based on found specifications
% --------------------------------------------------------------------------------
% R2=fWTFindPitch(WT ,4:2:25,Pref,1,Algo);
%WT.Spec.vSIMRef=R.vSIMRef;
%R2=fWTPowerCurve('BEM',WT ,4:0.5:25,1,Algo)
%figure, fplotCodesComparison('WS','CT',{R2},'','','',1,1,[],[],'','');

%%
%legds={}
%%
% disp('Watch out you are going to clear the workspace')
% pause
% load('WTSimRefWorkspace.mat')
% %%
% figure, fplotCodesComparison('r','CTloc',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),axis ij,box on,colorbar
% figure, fplotCodesComparison('r','a',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),axis ij,box on,colorbar
% figure, fplotCodesComparison('r','Gamma',R.PowerCurveData','','','',1,1,[],[],'','',R.WS,'WS'),box on,colorbar
%% Wind Farm control activities


end

