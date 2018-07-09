function [ dp ] = fdp( test_pitch )
%FDP Summary of this function goes here
%   Detailed explanation goes here

        Sim.Run.PITCH=test_pitch;
        [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)','BEM'));
        dp=Result.Power-Pref;
end

