function [Power,Pitch,RPM,WS,Ud,BEM_a,BEM_r,BEM_CP,BEM_CT] = fControllerRef( derating,U)
%FCONTROLLERREF Summary of this function goes here
%   Detailed explanation goes here
global PATH
PATH.BEM='/home/ewan/Work/PostDoc/WindfarmControl/dev/BEM_v05/BEM/';
PATH.WTlib='/home/ewan/Work/PostDoc/WindfarmControl/dev/BEM_v05/WTlib/';
PATH.Wind='/home/ewan/Work/PostDoc/WindfarmControl/dev/BEM_v05/Wind/';
PATH.DATA_WT='/home/ewan/Work/PostDoc/WindfarmControl/dev/BEM_v05/WT-data/';

require('BEM','v05');
require('WTlib','v06');
require('Wind','v01');

%%
% sWT='SB2'; Format='xblade';
% sWT='Riso10MW'; Format='hawc';
% sWT='NREL5MW'; Format='hawc';
% sWT='NTK500'; Format='hawc';
sWT='NY2'; Format='hawc';

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
% WT.Nacelle.tilt=5.;  % NREL 5MW
% WT.Rotor.cone=2.5; % NREL 5MW
% WT.Nacelle.tilt=5.;  % NY2
% WT.Rotor.cone=2.7; % NY2
WT.Nacelle.tilt=0.0;
WT.Rotor.cone=0.0;

%%
Opts.TipSpeedMax=(WT.Rotor.R*pi*15.947)/30; % Ny2

Opts.OmegaMin=7.6394*2*pi/60; % NY2
% Opts.TipSpeedMax=80; % NREL5MW
% Opts.OmegaMin=6.9*2*pi/60; % NREL5MW
% Opts.bOptimBelowPref=true;
Opts.bOptimBelowPref=false;

Opts.vWS_out=4:0.5:25;
R=fWTPowerCurve('BEM',WT ,Opts.vWS_out,0,Algo);

WScurve=R.WS;
Powcurve=R.Power;
RPMcurve=R.RPM;
PITCHcurve=R.PITCH;

%%

Pd=max(R.Power)*1000*derating;
RR=fWTPowerCurve('BEM',WT ,U,0,Algo);
RR.Power
Prefs = interp1(WScurve,Powcurve,U);
% disp(Prefs)
if Pd >= WT.Spec.P_rated
    warning('P derated above Prated')
else
    [~,iP]=min(abs(Powcurve-WT.Spec.P_rated*0.001));
    Powcurve=Powcurve(1:iP);
    WScurve=WScurve(1:iP);
    RPMcurve=RPMcurve(1:iP);
    PITCHcurve=PITCHcurve(1:iP);
%     Pd=min(Pd,Prefs*1000)
    Ud   = interp1(Powcurve,WScurve,Pd*0.001); % <<<<<<< POWER NOT UNIQUE
    RPMd = interp1(WScurve,RPMcurve,Ud);
    PITCHd = interp1(WScurve,PITCHcurve,Ud);
end

% Setting up the derated "WS-RPM" curve based on vSIMRef and RPM_derated%

% vWSd=linspace(3,25,100);
%evaluate at one speed
vWSd=U;
[~,iUd]=min(abs(vWSd-Ud));

vRPM=interp1(WScurve,RPMcurve,vWSd);
vRPM(iUd:end)=RPMd;

vPITCH=interp1(WScurve,PITCHcurve,vWSd);
vPITCH(iUd:end)=PITCHd;

vPITCHRef= zeros(1,length(vWSd)); % Unknown

vSIMRef=[vWSd(:) vRPM(:) vPITCHRef(:)];
% vSIMRef=[vWSd(:) vRPM(:) vPITCH];
% vSIMRef=[vWSd(:) vRPM(:) PITCHd];
%%
disp(vSIMRef)
disp(vWSd)
disp(Pd)

R3=fWTFindPitchEval(WT, vSIMRef, vWSd,Pd,0,Opts.bOptimBelowPref,Algo);

Power=R3.Power;
Pitch=R3.PITCH;
RPM=R3.RPM;
WS=R3.WS;
BEM_a=R3.PowerCurveData.a;
BEM_r=R3.PowerCurveData.r;
BEM_CP=R3.PowerCurveData.CP;
% BEM_CPloc=R3.PowerCurveData.CPloc;
BEM_CT=R3.PowerCurveData.CT;


end

