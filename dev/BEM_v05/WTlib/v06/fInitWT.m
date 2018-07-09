function [ WT ] = fInitWT( varargin )
% sWT, sFormat,sPath, (Opts) ,(Files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Documentation - OLD DOCUMENTATION FROM fINIT BEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function initialize as much parameters as possible given the input
% files. These parameters are different for one format to the other. So one
% should be aware of them together with the values set by default in
% InitDefault.

% Rotor.r :
% --------
% this is the radial position along the blade span starting from
% the rotational axis and going along the coned line. It has values between
% Rotor.rhub and rotor.R (included or excluded)
% if Opts.KeepOriginalR, then the r included in the param files is used
% if Opts.ForcedR then Opts.ForcedR is used for r
% else, a linear distribution of Algo.Ngrid points is used

% To use the AeroDyn integration method, r_hub and R should not be part of
% Rotor.r. They are always part of Rotor.rfull though.
% In the aerodyn integration method, Rotor.dr is used, which is the length
% of each blade elements

% Coning:
%---------
% If the rotor is coned, the Rotor.Rconed < Rotor.R
% This should be used for the tip-speed ratio abd for the swept area


WT=fInitWTDefault();

if(nargin>0)
    % for i=1:floor(nargin/2)
    %     switch(varagin(2*(i-1)+1))
    %     end
    % end
    %% Arguments
    sWT=varargin{1};
    sFormat=lower(varargin{2});
    sPath=varargin{3};
    if(nargin>3)
        Opts=varargin{4};
    else
        Opts=[];
    end
    if(nargin>4)
        Files=varargin{5};
        % Spec file
        Ispec =whichfile(Files, '(Spec\.dat)$');
        if(~isempty(Ispec))
%             disp(['Reading Spec File:',Files{Ispec}])
            fReadSpec(Files{Ispec});
        end
    end
    %% Common stuff
    WT.Name=sWT;
    WT.Sources.Format=sFormat;
    %% Load Aero Files
    % finds the wind turbine folder
    Aero =fFindSimFolder([sPath],sWT,0,0,0,'');
    if(~isfield(Aero,'WTPath'))
        warning('Note: no Aero Wind turbine Path')
    else
        WT.Sources.WTPath=Aero.WTPath;
        %%% XBlade
        switch sFormat
            case 'xblade'
                buff=dir([Aero.DataPath Aero.WTPath '*_param.txt']);
                Files{1}=[Aero.DataPath Aero.WTPath buff.name];
                buff=dir([Aero.DataPath Aero.WTPath '*.pc']);
                Files{2}=[Aero.DataPath Aero.WTPath buff.name];
                Opts.Extended=~isempty(regexpi(buff.name,'extended'));
                WT=fInitXblade(Files,WT,Opts);
            case 'hawc'
                if ~exist('Files')
                    disp( Aero.DataPath)
                    
                    buff=dir([Aero.DataPath Aero.WTPath '*.htc']);
                    Files{1}=[Aero.DataPath Aero.WTPath buff(1).name];
                    buff=dir([Aero.DataPath Aero.WTPath '*Spec.dat']);    
                    Files{2}=[Aero.DataPath Aero.WTPath buff(1).name];
                end
                Algo.Format='hawc';
                WT=fInitHawc(Files,WT,Opts);
            case 'flex'
%                 keyboard
                if ~exist('Files')
                    buff=dir([Aero.DataPath Aero.WTPath '*BladeGeometry.dat']);
                    Files{1}=[Aero.DataPath Aero.WTPath buff(1).name];
                    buff=dir([Aero.DataPath Aero.WTPath '*BladeProfiles.dat']);
                    Files{2}=[Aero.DataPath Aero.WTPath buff(1).name];
                    buff=dir([Aero.DataPath Aero.WTPath '*Spec.dat']);
                    Files{3}=[Aero.DataPath Aero.WTPath buff(1).name];
                end
                Algo.Format='flex';
                WT=fInitFlex(Files,WT,Opts);
            case 'wtperf'
                Algo.Format='wtperf';
                WT=fInitWTPerf(Files,Opts);
            otherwise
                warning('Unknown format for rotor specification initialization')
        end
    end
end



WT.Spec.vSIM=[];

%% Specific stuff
% Geometry
if(isequal(sWT,'SB2'))
    WT.Spec.P_rated= 2000*1000; 
    WT.Spec.vSIM=[
        06 10.7 -2
        8 14.3 -2
        10 15.9 -1.4
        12 16 4.5
        15 16 10.4
%         20 16 17.6 
];


WT.Spec.vSIMRef=[
   3.000000     5.500000    -2.000000
    4.000000     7.250000   -2.000000
    5.000000     8.250000   -2.000000
    6.000000    11.250000   -2.000000
    7.000000    13.000000   -2.000000
    8.000000    14.750000   -2.000000
    9.000000    16.000000   -2.000000
   10.000000    16.000000   -1.500000
   11.000000    16.000000    1.000000
   12.000000    16.000000    4.000000
   13.000000    16.000000    6.500000
   14.000000    16.000000    8.500000
   15.000000    16.000000   10.400000
   20.000000    16.000000   17.600000   ];


    WT.DefaultEnvironment.rho=1;
    WT.DefaultEnvironment.KinVisc=18.3*10^-6;
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'B49clean'))
    WT.Spec.vSIM=[
        06 10.7 -2
        8 14.3 -2
        10 15.9 -1.4
        12 16 4.5
        15 16 10.4
%         20 16 17.6 
        ];
    WT.DefaultEnvironment.rho=1;
    WT.DefaultEnvironment.KinVisc=18.3*10^-6;
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'Manu'))
    WT.Spec.vSIM=[06 4.2 0];
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'Mexico'))
    WT.Spec.vSIM=[10 424.4 0];
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'NTK500'))
    WT.Spec.P_rated= 500*1000; 
    WT.Spec.vSIM=[06 27.1 0];
    WT.Spec.vSIMRef=[
  5.0 27.1   0.00 
  6.0 27.1   0.00 
  7.0 27.1   0.00 
  8.0 27.1   0.00 
  9.0 27.1   0.00 
 10.0 27.1   0.00 
 11.0 27.1   0.00 
 12.0 27.1   0.00 
 13.0 27.1   0.00 
 14.0 27.1   0.00 
 15.0 27.1   0.00 
 16.0 27.1   0.00 
 17.0 27.1   0.00 
 18.0 27.1   0.00 
 19.0 27.1   0.00 
 20.0 27.1   0.00 
 21.0 27.1   0.00 
 22.0 27.1   0.00 
 23.0 27.1   0.00 
 24.0 27.1   0.00 
 25.0 27.1   0.00 ];

    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'NTK500p'))
    WT.Spec.P_rated= 500*1000; 
    WT.Spec.vSIM=[06 27.1 0];
    WT.Spec.vSIMRef=[
  5.0 27.1   0.00 
  6.0 27.1   0.00 
  7.0 27.1   0.00 
  8.0 27.1   0.00 
  9.0 27.1   0.00 
 10.0 27.1   0.00 
 11.0 27.1   0.00 
 12.0 27.1   0.00 
 13.0 27.1   2.19 
 14.0 27.1   6.66 
 15.0 27.1   9.72 
 16.0 27.1  12.27 
 17.0 27.1  14.51 
 18.0 27.1  16.56 
 19.0 27.1  18.45 
 20.0 27.1  20.23 
 21.0 27.1  21.91 
 22.0 27.1  23.50 
 23.0 27.1  25.01 
 24.0 27.1  26.46 
 25.0 27.1  27.84 ];

    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'NRELShen'))
    WT.Spec.vSIM=[
        07 71.93 5
        10 71.93 5
        13 71.93 5
        15 71.93 5        ];
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'V90'))
    WT.Spec.P_rated= 3000*1000; 
    WT.Spec.vSIM=[ ];
    WT.Spec.iSimDef=1;
    WT.Spec.vSIMRef=[];
elseif(isequal(sWT,'NREL5MW_default'))
    WT.Spec.P_rated= 5000*1000; 
    WT.Spec.vSIM=[ 7.0000    8.4883         0];
    WT.Spec.iSimDef=1;
    WT.Spec.vSIMRef=[
    4.0000    4.8504         0
    5.0000    6.0630         0
    6.0000    7.2757         0
    7.0000    8.4883         0
    8.0000    9.7009         0
    9.0000   10.9135         0
   10.0000   10.9135         0
   11.0000   10.9135         0
   12.0000   10.9135    4.3225
   13.0000   10.9135    7.1709
   14.0000   10.9135    9.4192
   15.0000   10.9135   11.3683
   16.0000   10.9135   13.1327
   17.0000   10.9135   14.7648
   18.0000   10.9135   16.2933
   19.0000   10.9135   17.7314
   20.0000   10.9135   19.1029
   21.0000   10.9135   20.4281
   22.0000   10.9135   21.7058
   23.0000   10.9135   22.9328
   24.0000   10.9135   24.1171
   25.0000   10.9135   25.2636];
elseif(isequal(sWT,'NREL5MW'))
    WT.Spec.P_rated= 5000*1000; 
    WT.Spec.vSIM=[ 7.0000    8.4883         0];
    WT.Spec.iSimDef=1;
    WT.Spec.vSIMRef=[
    3.0000    6.9000    1.5843
    4.0000    6.9000    1.2381
    5.0000    6.9000    0.2676
    6.0000    7.2757   -0.5404
    7.0000    8.4883   -0.5404
    8.0000    9.7009   -0.5404
    9.0000   10.9135   -0.5404
   10.0000   12.1261   -0.5404
   11.0000   12.1261   -0.7100
   12.0000   12.1261    4.6275
   13.0000   12.1261    7.0817
   14.0000   12.1261    9.0442
   15.0000   12.1261   10.7616
   16.0000   12.1261   12.3242
   17.0000   12.1261   13.7751
   18.0000   12.1261   15.1359
   19.0000   12.1261   16.4230
   20.0000   12.1261   17.6617
   21.0000   12.1261   18.8653
   22.0000   12.1261   20.0292
   23.0000   12.1261   21.1496
   24.0000   12.1261   22.2335
   25.0000   12.1261   23.2875];
elseif(isequal(sWT,'NY2'))
    WT.Spec.P_rated= 2300*1000; 
    WT.Spec.vSIM=[7.00     11.85961028  -1.5    ];
    WT.Spec.iSimDef=1;
    WT.Spec.vSIMRef=[4.00     6.77692016  -1.5       
4.25     7.20047767  -1.5       
4.50     7.62403518  -1.5       
4.75     8.04759269  -1.5       
5.00     8.4711502   -1.5       
5.25     8.89470771  -1.5       
5.50     9.31826522  -1.5       
5.75     9.74182273  -1.5       
6.00     10.16538024  -1.5       
6.25     10.58893775  -1.5       
6.50     11.01249526  -1.5       
6.75     11.43605277  -1.5       
7.00     11.85961028  -1.5       
7.25     12.28316779  -1.5       
7.50     12.7067253   -1.5       
7.75     13.13028281  -1.5       
8.00     13.55384031  -1.5       
8.25     13.97739782  -1.5       
8.50     14.40095533  -1.5       
8.75     14.82451284  -1.5       
9.00     15.24807035  -1.5       
9.25     15.67162786  -1.5       
9.50     15.947       -1.5       
9.75     15.947       -1.5       
10.00     15.947       -1.5       
10.25     15.947       -1.5       
10.50     15.947       -1.5       
10.75     15.947       -1.5       
11.00     15.947        1.20448849
11.25     15.947        2.47250561
11.50     15.947        3.42626409
11.75     15.947        4.25016791
12.00     15.947        4.99108715
12.25     15.947        5.66961967
12.50     15.947        6.29923152
12.75     15.947        6.89076581
13.00     15.947        7.45191567
13.25     15.947        7.98823228
13.50     15.947        8.50364739
13.75     15.947        9.00097951
14.00     15.947        9.48220989
14.25     15.947        9.94911337
14.50     15.947       10.40290607
14.75     15.947       10.84508107
15.00     15.947       11.27648448
15.25     15.947       11.69808982
15.50     15.947       12.11066831
15.75     15.947       12.51511524
16.00     15.947       12.91195496
16.25     15.947       13.30163338
16.50     15.947       13.6847893 
16.75     15.947       14.06181858
17.00     15.947       14.4332778 
17.25     15.947       14.79949631
17.50     15.947       15.16041755
17.75     15.947       15.5165319 
18.00     15.947       15.86793277
18.25     15.947       16.21500877
18.50     15.947       16.55606561
18.75     15.947       16.893017  
19.00     15.947       17.22621927
19.25     15.947       17.55614276
19.50     15.947       17.88265862
19.75     15.947       18.20580651
20.00     15.947       18.52604004
20.25     15.947       18.84318651
20.50     15.947       19.15719251
20.75     15.947       19.46857929
21.00     15.947       19.77719174
21.25     15.947       20.08297159
21.50     15.947       20.38605063
21.75     15.947       20.68670497
22.00     15.947       20.98503649
22.25     15.947       21.28094783
22.50     15.947       21.57444652
22.75     15.947       21.8656929 
23.00     15.947       22.15477889
23.25     15.947       22.44179025
23.50     15.947       22.72630348
23.75     15.947       23.00883895
24.00     15.947       23.28906478
24.25     15.947       23.56726957
24.50     15.947       23.84347192
24.75     15.947       24.11767314
25.00     15.947       24.38993115];
elseif(isequal(sWT,'NREL'))
    WT.Spec.vSIM=[
        07 71.93 5
        10 71.93 5
        13 71.93 5
        15 71.93 5        ];
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'NREL_extended'))
    WT.Spec.vSIM=[06 72.1 5];
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'Riso10MW'))
    WT.Spec.P_rated= 11000*1000; 
%     WT.Spec.vSIM=[09 7.72 0];
    WT.Spec.vSIM=[10 8.06 0]; % Frederik AIAA paper 
    WT.Spec.iSimDef=6;
    WT.Spec.vSIMRef=[ % Frederik AIAA paper
        5.0  6.000  1.966
        6.0  6.000  0.896
        7.0  6.000  0.000
        8.0  6.426  0.000
        9.0  7.229  0.000
        10.0 8.060  0.000 % adapated RPM a bit
        11.0 8.836  0.000
        12.0 9.600  4.502
        16.0 9.600  12.499
        20.0 9.600  17.618
        25.0 9.600  22.975 ];
%     WT.Spec.vSIMRef=[
%     4.0000    3.4298         0
%     5.0000    4.2873         0
%     6.0000    5.1447         0
%     7.0000    6.0022         0
%     8.0000    6.8597         0
%     9.0000    7.7171         0
%    10.0000    7.7171         0
%    11.0000    7.7171         0
%    12.0000    7.7171    4.3166
%    13.0000    7.7171    7.1653
%    14.0000    7.7171    9.4136
%    15.0000    7.7171   11.3628
%    16.0000    7.7171   13.1272
%    17.0000    7.7171   14.7594
%    18.0000    7.7171   16.2879
%    19.0000    7.7171   17.7262
%    20.0000    7.7171   19.0977
%    21.0000    7.7171   20.4229
%    22.0000    7.7171   21.7006
%    23.0000    7.7171   22.9276
%    24.0000    7.7171   24.1119
%    25.0000    7.7171   25.2585 ];
% 


elseif(isequal(sWT,'SB1'))
    WT.Spec.vSIM=[06 09 -1.3;  10 8 -2.7;10 12 -2.7; 12 12 6.85;14 12 10.85;  16 12 15.7];
    WT.Spec.iSimDef=3;
elseif(isequal(sWT,'SB1_tip')) %new tip
    WT.Spec.vSIM=[06 09 -2.5;  06 09 -2.8 ; 06 09 -3.25 ; 10 12 -2.7];
    WT.Spec.iSimDef=4;
elseif(isequal(sWT,'SB1_swept'))
    %WT.Spec.vSIM=[10 06 -2.7; 10 08 -2.7; 10 12 -2.7; 10 14 -2.7];
    WT.Spec.vSIM=[10 12 -2.7; 10 14 -2.7];
    WT.Spec.iSimDef=1;
elseif(isequal(sWT,'Tjaere'))
    WT.Spec.vSIM=[8.7 22 0];
    WT.Spec.iSimDef=1;
    WT.Sepc.P_rated=2200*1000;
elseif(isequal(sWT,'FlatTUD'))
    WT.Spec.vSIM=[05.5 700.2817 -3.1];
    WT.Spec.iSimDef=1;
    WT.Rotor.nB=2;
elseif(isequal(sWT,'elliptic'))
    WT.Spec.vSIM=[1 0 0];
    WT.Rotor.nB=1;
    WT.Rotor.R=10;
elseif(isequal(sWT,'simpleblade'))
    WT.Spec.vSIM=[10 0 0];
    WT.Rotor.nB=3;
    WT.Rotor.R=10;
elseif(isequal(sWT,'flatplate'))
    WT.Spec.vSIM=[10 0 0];
    WT.Rotor.nB=1;
    WT.Rotor.R=10;
elseif(isequal(sWT,'trapezoid'))
    WT.Spec.vSIM=[10 0 0];
    WT.Rotor.nB=1;
    WT.Rotor.R=2;
elseif(isequal(sWT,'B49XLL'))
    WT.Spec.vSIM=[
        1.00000000e+000   1.46880707e+000 0
        3.00000000e+000   4.40642120e+000 0
        5.00000000e+000   7.34403534e+000 0
        7.00000000e+000   1.02816495e+001 0
        9.00000000e+000   1.32192636e+001 0
        1.10000000e+001   1.61568777e+001 0
        1.30000000e+001   1.86000000e+001 0
        1.50000000e+001   1.86000000e+001 0
        1.70000000e+001   1.86000000e+001 0
        1.90000000e+001   1.86000000e+001 0
        2.10000000e+001   1.86000000e+001 0
        2.30000000e+001   1.86000000e+001 0
        2.50000000e+001   1.86000000e+001 0 ];
    WT.Spec.iSimDef=4;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(WT.Rotor.cone~=0)
    warning('This rotors has some coning!')
end
if(WT.Rotor.R==-1)
    warning('The rotor length has not been initialized - I set it to the max!')
    WT.Rotor.R=max(WT.Sources.Rotor.r);
end
if(WT.Rotor.nB==-1)
    warning('The rotor length has no blades... I set it to 3!')
    WT.Rotor.nB=3;
end


if(sum(isnan(WT.Sources.Rotor.twist))>0)
    warning('WT.Sources.Rotor.Twist has NaN values!')
    kbd
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Element used for the estimate of khi (at r/R=70%)
WT.Rotor.R_coned =WT.Rotor.R*cosd(WT.Rotor.cone);   %  Rotor coned radius [m]
WT.Rotor.SweptArea=pi*(WT.Rotor.R*cosd(WT.Rotor.cone))^2;
WT.Spec.TSR_rated=WT.Spec.Omega_rated*WT.Rotor.R*cosd(WT.Rotor.cone)/WT.Spec.V_rated;

end

