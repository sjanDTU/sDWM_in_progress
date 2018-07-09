function [ mp ] = fmin_minus_p( test_pitch)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
       Sim.Run.PITCH=test_pitch;
       [ Result ] = eval(sprintf('fRun%s(WT,Sim,Wind,Algo)','BEM'));
       mp=-Result.Power;

end

