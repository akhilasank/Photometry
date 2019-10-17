% function params = HL_FP_Wheel_default_params
% generate all the default parameters for Wheel speed encoder processing
%
% Haixin Liu 2019-9
% taken from Pratiks processParams.m script
function params = HL_FP_Wheel_default_params


params.beh.radius = 9.8; %Radius of the wheel used. Note it can be meters or centimeters. Just keep track of your units
params.beh.winSize = 0.5; %This is the window size for the moving median filter applied to unwrapped encoder data 500ms windows work well
%Onset/Offset Parameters
%Movement Onset and Offset Parameters
params.beh.velThres = 4; %(same units as radius)/s
params.beh.minRunTime = 4; %Threshold for minimum time spent running for movement bouts (in seconds)
params.beh.minRestTime = 4; %Threshold for minimum time spent rest for movement bout (in seconds)
params.beh.finalOnset = 0; %Boolean value -- Decides if you want to include or exlcude the final 
% onset if the acquisition ends before the offset
params.beh.timeThres = 4; %Make sure a bout is above a certain time-length
params.beh.timeBefore = 4; %Time to display preceding movement onset and offset
params.beh.timeAfter = 4; %Time to display following movement onset and offset

%Rest Onset and Offset Parameters
params.beh.minRestTime_rest = 4;
params.beh.minRunTime_rest = 1;
params.beh.velThres_rest = 2;
params.beh.timeThres_rest = 4;
params.beh.timeShift_rest = 0.5;

