% function data = HL_FP_process_Wheel(data, ts, rawFs, params)
% to process running wheel data genrated by speed encoder
% INPUT: 
%       data: vector, nx1, of the raw voltage recording
%       params: -optional, -struct: params.beh.radius 
%                                          beh.velThres 
%                                          beh.winSize
%                                          beh.finalOnset
% OUTPUT:
%       data: 
%
% Haixin Liu 2019-9
% revised from processBeh.m in Process Function folder
% NOT WORKING YET
function data = HL_FP_process_Wheel(data,params)

error('Function not in use for now')
if nargin < 2 % default
    disp('Using default wheel parameters:');
    
else
radius = params.beh.radius; velThres = params.beh.velThres;
winSize = params.beh.winSize;
finalOnset = params.beh.finalOnset;
end


nAcq = length(data.acq);
for n = 1:nAcq
    if (isfield(data.final(n),'wheel'))
        if (~isempty(data.final(n).wheel))
            wheel = data.final(n).wheel;
        else
            wheel = data.acq(n).wheel;
        end
    else
        wheel = data.acq(n).wheel;
    end
    Fs = data.acq(n).Fs;
    if params.dsRate ~= 0
       wheel = downsample(wheel,params.dsRate);
       Fs = Fs/params.dsRate;
    end
    data.final(n).wheel = wheel;
    data.final(n).Fs = Fs;
    vel = getVel(wheel,radius,Fs,winSize);
    minRest = params.beh.minRestTime * Fs; minRun = params.beh.minRunTime * Fs;
    [onsets,offsets] = getOnsetOffset(abs(vel),velThres,minRest,minRun,finalOnset);
    data.final(n).vel = vel;
    data.final(n).beh.onsets = onsets;
    data.final(n).beh.offsets = offsets;
end
