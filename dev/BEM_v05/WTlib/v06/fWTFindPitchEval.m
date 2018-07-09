function [R]= fWTFindPitchEval(WT,vSIMRef,vWS,Pref,bPlot,bOptimBelowPref,Algo)
Code='BEM';

% if(~isfield(WT.Spec,'vSIMRef'))
%     if(~isfield(WT.Spec,'vSIM'))
%         error('Provide either WT.Spec.vSIM or vSIMRef')
%     else
%         vSIMRef=WT.Spec.vSIM;
%     end
% end
% if min(vSIMRef(:,1))>min(vWS) || max(vSIMRef(:,1))<max(vWS)
%     warning('Not enough input data in SimRef, will have to extrapolate');
% end
% if(length(unique(vSIMRef(:,2)))==1)
%     fprintf('!!! Rpm is constant in specs\n');
% end
% pitch = interp1(vSIMRef(:,1), vSIMRef(:,3),vWS,'cubic','extrap'); % changed from cubic to linear for debugging
% rpm   = interp1(vSIMRef(:,1), vSIMRef(:,2),vWS,'cubic','extrap');
% if(length(unique(rpm))==1)
%     fprintf('!!! Using constant rpm\n');
% end
% vSIM=[vWS(:) rpm(:) pitch(:)];
vSIM=vSIMRef;
[ Sim ]  = fInitSim(WT);
[ Wind ] = fInitWind(  ); 
% keyboard

optim_opts = optimset('TolX',5e-4);
% optim_opts = optimset('TolX',5e-2,'Display','final');

pitch=zeros(1,length(vWS));
power=zeros(1,length(vWS));
for i=1:length(vWS)
    pitch(i)=vSIM(i,3);
    Sim.Run=fSetRun(vSIM(i,:));
    Wind=fSetWind(Wind,Sim);
    %first run for fun
    [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
    P_before=Result.Power;
    DP = 0;
    if Result.Power ~= Pref
        if i>1 
            % we use the previous pitch+6.1 deg as a starting point for fzero
            pitch(i) = fzero(@fdp,pitch(i-1)+6.1,optim_opts);
        else
            disp('yes')
            pitch(i) = fzero(@fdp,1.1,optim_opts);
        end
        Sim.Run.PITCH=pitch(i);
        disp('chosen pitch is')
        disp(pitch(i))
        [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
    else
        if bOptimBelowPref
            % Optimization of power below P rated (especially if Omega Min)
            pitch(i) = fminbnd(@fmin_minus_p,pitch(i)-5,pitch(i)+5,optim_opts); 
            Sim.Run.PITCH=pitch(i);
            [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
            DP       = (Result.Power-P_before)/1000;
        end
    end
    R.PowerCurveData(i)=Result;
    power(i) = Result.Power/1000 ;   
    fprintf('U = %5.1f  - P = %7.1f kW - pitch = %5.2f deg - Rpm = %4.1f (DP = %7.1f kW)\n',vWS(i),Result.Power/1000,pitch(i),Sim.Run.RPM,DP);
end
R.PITCH=pitch;
R.Power=Result.Power;
R.WS=vWS;
R.RPM=vSIM(:,2);
R.SIMRef=[R.WS(:) R.RPM(:) R.PITCH(:) R.PITCH(:)*0];

% nested function
function dp = fdp(test_pitch)
        Sim.Run.PITCH=test_pitch;
%         Sim.Run.PITCH=max(test_pitch,-2.5);
%         Sim.Run.PITCH=min(test_pitch,35);
        [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
        dp=Result.Power-Pref;
        disp([Sim.Run.PITCH,Result.Power/1000,Pref/1000,dp])
end
function mp = fmin_minus_p(test_pitch)
        Sim.Run.PITCH=test_pitch;
        [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)',Code));
        mp=-Result.Power;
end


if(bPlot)
    Codes={R}; legds={Code};
    colrs=fColrs(1:4);
    sty={'-','+-','--'};
    figure, fplotCodesComparison('WS','PITCH',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','RPM',Codes,legds,colrs,sty,1,1,[],[],'','')
    figure, fplotCodesComparison('WS','Power',Codes,legds,colrs,sty,1,1,[],[],'','')
    dispatchFigs(1);
end

end % end function
