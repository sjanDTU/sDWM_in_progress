function [RES] = fBEMsteady(WT,Sim,Wind,Algo) %#codegen
cone=WT.Rotor.cone;

V0=Wind.V0(3)*cosd(cone);
VHubHeight=Wind.V0(3);

nB=WT.Rotor.nB;
ne=WT.Rotor.ne;
r=WT.Rotor.r;
dr=WT.Rotor.dr;
rhub=WT.Rotor.rhub;
R=WT.Rotor.R;

chord=WT.Rotor.chord;
twist=WT.Rotor.twist;
cone=WT.Rotor.cone;
Rotor=WT.Rotor;

sigma=chord*nB./(2*pi*r*cosd(cone));  %!!! CONE

% Environment
rho=Sim.rho;
KinVisc=Sim.KinVisc;

pitch=Sim.Run.PITCH;
omega=Sim.Run.RPM*2*pi/60;
lambda_r=omega*r*cosd(cone)/V0;       %!!! cone
lambda=omega*R*cosd(cone)/V0;         %!!! cone


%algorithm internal paramter
nbIt=Algo.nbIt;
aTol=Algo.aTol;


% initialize vectors
Pn=zeros(1,ne);
Pt=zeros(1,ne);
RES.Vrel=zeros(1,ne);
RES.Re=zeros(1,ne);
RES.F=zeros(1,ne);
RES.Fshen=zeros(1,ne);
RES.a=zeros(1,ne);
RES.aprime=zeros(1,ne);
RES.phi=zeros(1,ne);
RES.alpha=zeros(1,ne);
RES.Cl=zeros(1,ne);
RES.Cd=zeros(1,ne);
RES.Cn=zeros(1,ne);
RES.Ct=zeros(1,ne);
RES.CT=zeros(1,ne);
RES.nIt=zeros(1,ne);


if(isequal(Algo.BEM.TipLossMethod,'GoldsteinSimple'))
    l_bar=1/lambda;
    vx=linspace(0,1,3*length(r));
    FGoldstein=fTipLossGoldsteinOkulov( l_bar,nB,vx );
    FGoldstein(FGoldstein<0)=0;
    FGoldstein(isnan(FGoldstein))=1;    
    %     figure
    %     plot(vx,FGoldstein)

end


% for e=ne:-1:1 % !!!!!!!!!!!!!!!!!! backward for Xu and Sankar correction
for e=1:ne
    %initialization
    a=0.3*0;
    aprime=0.01*0;
    Ftip_previous=1;
    %RES loop
    for i=1:nbIt
        % --------------------------------------------------------------------------------
        % --- Step 0: Relative wind
        % --------------------------------------------------------------------------------

        % --------------------------------------------------------------------------------
        % --- Step 1: Wind Components
        % --------------------------------------------------------------------------------
        Ut = omega*r(e)*(1+aprime);
        Un = V0*(1-a);
        Vrel_norm = sqrt(Un.^2+Ut.^2); %more stable than  Vrel=V0*(1-a)/sind(phi);
        Re=Vrel_norm*chord(e)/KinVisc/10^6; % Reynolds number in Millions
%         disp([num2str(Ut),' ',num2str(Un),' ',num2str(Vrel_norm),' ',num2str(Re)]) 
        if(imag(Un)+imag(Ut)~=0)
%             warning('crash');
%             keyboard
        end
        % --------------------------------------------------------------------------------
        % --- Step 2: Flow Angle
        % --------------------------------------------------------------------------------
        phi = atan2(Un,Ut)*180/pi; %(this returns a positive flow angle) [deg]
        if(imag(phi)~=0)
            disp(['Algorithm did not converge :',sprintf('lambda=%.2f beta=%.2f V0=%.2f r=%.2f - it %d',lambda,twist(e),V0,r(e),i)])
            phi=0;
            break
        end


        % --------------------------------------------------------------------------------
        % --- Tip loss
        % --------------------------------------------------------------------------------
        Ftip=a*0+1;
        Fperf=a*0+1; 
        if(Algo.BEM.bTipLoss)
            % tip loss correction
            switch(Algo.BEM.TipLossMethod)
                case 'Glauert'
                    if(sind(phi)>0.01   ) %singularity if sin phi close to zero and phi negative....
                        Ftip=2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi))));
                    end
                case 'Prandtl'
                    Ftip=2/pi*acos(exp(-nB/2*(1-lambda_r(e)/lambda) *sqrt(1+lambda^2 )));
                case 'GoldsteinSimple'
                    Ftip=interp1(vx,FGoldstein,r(e));
                    Ftip=max(min(Ftip,1),0);
                case 'XuSankar'
                    if(abs(sind(phi))>0.01) %singularity if sin phi close to zero
                        if(r(e)/R>=0.7)
                            FGl=2/pi*acos(exp(-nB/2*(R-r(e))/(R*sind(phi))));
                            Ftip=0.5*(FGl^0.85+0.5);
                        else
                            i07=whichvalue(r,0.7*R);
                            FGl=2/pi*acos(exp(-nB/2*(R-0.7*R)/(R*sind(RES.phi(i07)))));
                            F=0.5*(FGl^0.85+0.5);
                            Ftip=1-r(e)/R*(1-F)/0.7 ;
                        end
                    end
                case 'Lindenburg'
                    Ftip=2/pi*acos(exp(-nB/2*(1-lambda_r(e)/lambda)*sqrt(1+lambda_r(e)^2*   ( (1+2*sqrt(Ftip_previous)*aprime/2)/(1-sqrt(Ftip_previous)*a/2) )^2    )));
                case 'Shen'
                    if(abs(sind(phi))>0.01) %singularity if sin phi close to zero
                        Ftip=2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi))));
                        g=exp(-0.125*(nB*lambda-21))+0.1;
                        Fperf  =2/pi*acos(exp(-nB/2*(R-r(e))/(r(e)*sind(phi)) *g));
                    end
                case 'PrescribedWake'
                    error('I should not be run by this script but by the dedicated standalone one')                                      
                case 'TipLossDB'
                    Fperf=2/pi*acos(exp(-63*(1-r(e)/R)));
                otherwise
                    error('tip loss method unknown')
            end
        end
        % --------------------------------------------------------------------------------
        % --- Hub loss
        % --------------------------------------------------------------------------------
        Fhub=a*0+1;
        if(Algo.BEM.bHubLoss)
            %prandtl hub loss correction
            Fhub=2/pi*acos(exp(-nB/2*(r(e)-rhub)/(rhub*sind(phi))));
        end
        F=Ftip.*Fhub;
        
        % --------------------------------------------------------------------------------
        % --- Step 3: Angle of attack
        % --------------------------------------------------------------------------------
        if(Algo.bAlphaHacked)
            alpha=interp1(Rotor.AlphaHacked(:,1),Rotor.AlphaHacked(:,2),r(e));
            phi=alpha+(twist(e)+pitch);
        else
             alpha=phi-(twist(e)+pitch); %[deg]
        end
        
%         keyboard
        % --------------------------------------------------------------------------------
        % --- Step 4: Profile Data
        % --------------------------------------------------------------------------------
        if(sum(isnan(alpha))>0)
%             warning('Alpha bad')
%             keyboard
        end
%         disp([num2str(phi),' ',num2str(alpha),' ',num2str(twist(e)),' ',num2str(pitch)])

        [Cl Cd Cn Ct CnForAI CtForTI ] = fAeroCoeffWrap(e,alpha,phi,chord,Vrel_norm,Re,Fperf,WT,Algo); %Includes prescibed circulation
%         disp([num2str(Cl),' ',num2str(Cd),' ',num2str(Cn),' ',num2str(Ct),' ',num2str(CnForAI),' ',num2str(CtForTI)])
%         disp([num2str(Cl),' ',num2str(Cd)])
        % --------------------------------------------------------------------------------
        % --- Induction Factor in coordinate system
        % --------------------------------------------------------------------------------

        % --------------------------------------------------------------------------------
        % --- Step 5: Induction Coefficients
        % --------------------------------------------------------------------------------
        %a=(V0-Un)/V0;   %induction factor %abs for compatibility with unsteady code
        % Storing last values
        a_last=a;
        aprime_last=aprime;
        [ a aprime CT] = fInductionCoefficients(a_last,[0;0;Vrel_norm],Un,Ut,[0;0;V0],[0;0;V0],[0;0;-a_last*V0],omega,chord(e),F,Ftip,CnForAI,CtForTI,lambda_r(e),sigma(e),phi,Algo) ;

        % --------------------------------------------------------------------------------
        % --- Convergence Criteria
        % --------------------------------------------------------------------------------
        if (i>3 && abs(a-a_last)+abs(aprime-aprime_last)<aTol) %used to be on alpha
            break;
        end

        Ftip_previous=Ftip;

     
    end %end iterative loop for one element
%     disp([num2str(Cl),' ',num2str(Cd),' ',num2str(Cn),' ',num2str(Ct),' ',num2str(CnForAI),' ',num2str(CtForTI)])
    RES.nIt(e)=i;
    if(i==nbIt)
        fprintf('Maximum iterations reached : lambda=%.2f beta=%.2f V0=%.2f r=%.2f\n',lambda,twist(e),V0,r(e));
    end


    % --------------------------------------------------------------------------------
    % --- Step 6: Aerodynamic Forces PER LENGTH
    % --------------------------------------------------------------------------------
    %     L=0.5*rho*Vrel_norm.^2*chord(e)*Cl;
    %     D=0.5*rho*Vrel_norm.^2*chord(e)*Cd;
    %     Pn(e) = L*cosd(phi) + D*sind(phi);   %load normal to the rotor plane
    %     Pt(e) = L*sind(phi) - D*cosd(phi);   %load tangential to the rotor plane
    Pn(e) =   0.5*rho*Vrel_norm.^2*chord(e).*Cn;
    Pt(e) =   0.5*rho*Vrel_norm.^2*chord(e).*Ct;

    RES.Vrel(e)=Vrel_norm;
%     RES.Un(e)=Un;
%     RES.Ut(e)=Ut;
    RES.Re(e)=Re;
    RES.F(e)=F;
    RES.Fperf(e)=Fperf;
    RES.a(e)=a;
    RES.a_last(e)=a_last;
    RES.aprime(e)=aprime;
    RES.aprime_last(e)=aprime_last;
    RES.phi(e)=phi;
    RES.alpha(e)=alpha;
    RES.Cl(e)=Cl;% USED TO BE MULTIPLIED BY Fshen, why????
    RES.Cd(e)=Cd;% SAME 
    RES.Cn(e)=Cn;
    RES.Ct(e)=Ct;
    RES.CT(e)=CT;
end %loop on elements
%  keyboard
% Hawc convergence criteria norm(RES.aprime-RES.aprime_last)+norm(RES.a-RES.a_last)

RES.Pn=Pn;
RES.Pt=Pt;
RES.ThrLoc=dr.*Pn*cosd(cone);
RES.ThrLocLn=Pn*cosd(cone);
RES.CTloc=nB*RES.ThrLoc./(0.5*rho*VHubHeight^2*(2*pi*r.*cosd(cone).*dr)) ;

RES.TqLoc=dr.*r.*Pt*cosd(cone);
RES.TqLocLn=r.*Pt*cosd(cone);
RES.CQloc=nB.*RES.TqLoc./(0.5*rho*VHubHeight^2*(2*pi*r.*cosd(cone).*dr.*r*cosd(cone))) ;

%%%% Returning Aerodynamic Forces
if(isequal(WT.Sources.Format,'wtperf'))
    RES.Thrust = nB*sum(Rotor.dr.*(Pn*cosd(cone)));    %Rotor shaft thrust at t in Newton
    RES.Torque = nB*sum(Rotor.dr.*Pt.*(r*cosd(cone))); %Rotor shaft torque at t in Newton
else
    RES.Torque = nB*getTorqueFromBlade(r,Pt*cosd(cone),R);   %Rotor shaft torque at t in Newton
    RES.Thrust = nB*getThrustFromBlade(r,Pn*cosd(cone),R);   %Rotor shaft thrust at t in Newton
end
RES.Flap   = sum(dr.*(Pn*cosd(cone)).*(r-rhub));% approximate
RES.Edge   = sum(dr.*Pt.*(r*cosd(cone)).*(r-rhub));% approximate

RES.Power=omega*RES.Torque;
RES.CP=RES.Power/(0.5*rho*V0^3*WT.Rotor.SweptArea);
RES.CT=RES.Thrust/(0.5*rho*V0^2*pi*R^2);
RES.CQ=RES.Torque/(0.5*rho*V0^2*pi*R^3);
RES.Gamma=0.5*RES.Re.*RES.Cl*KinVisc*10^6;
RES.r=r;


RES.uia=V0*RES.a;
RES.uit=omega*r.*RES.aprime;

end

