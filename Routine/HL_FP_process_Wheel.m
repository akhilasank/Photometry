% function data = HL_FP_process_Wheel(data, ts, rawFs, params)
% to process running wheel data genrated by speed encoder
% INPUT: 
%       data: vector, nx1, of the raw voltage recording
%       Fs: sampling rate
%       params: -optional, -struct: params.beh.radius 
%                                          beh.velThres 
%                                          beh.winSize
%                                          beh.finalOnset
% OUTPUT:
%       vel:
%       wheel: location 
%       onsets: the indices of running epoch onset according to defined
%               params
%       offsets: the indices of offset
%
% Haixin Liu 2019-9
% revised from processBeh.m in Process Function folder

function [vel, wheel, distance, onsets,offsets, dsRate] = HL_FP_process_Wheel(data, Fs, dsRate, params)

% error('Function not in use for now')
if nargin < 4 % default
    disp('Using default wheel parameters:');
    params = HL_FP_Wheel_default_params;
    radius = params.beh.radius;
    velThres = params.beh.velThres;
    winSize = params.beh.winSize;
    finalOnset = params.beh.finalOnset;
    circum = 2*pi*params.beh.radius; 

else
    radius = params.beh.radius; 
    velThres = params.beh.velThres;
    winSize = params.beh.winSize;
    finalOnset = params.beh.finalOnset;
    circum = 2*pi*params.beh.radius; 
end

% process the rotatry encoder recording
% normalize the wheel data before processing the rotary encoder data
% empirical LEFT rig max = 4.9530 min 0.0178
temp_data = data; n_delete = 0;
if isempty(find(data>4))  % only pad when there is no min or max
    temp_data = [data; 4.9530];
    n_delete = n_delete+1;
end
if isempty(find(data<1))
    temp_data = [data; 0.0178];
    n_delete = n_delete+1;
end
temp_data = temp_data - min(temp_data);
temp_data = temp_data/max(temp_data);

% need to downsample first due to buggy performance with unwrap for high
% sampling or noisy data
% keep it the same as TLab, filter and downsample the data

% lpFilt = designfilt('lowpassiir','SampleRate',Fs,'FilterOrder',10,'HalfPowerFrequency',300);
% temp_data = filtfilt(lpFilt,temp_data);

if Fs > 1000
% filt = designfilt('lowpassiir','SampleRate',Fs,'HalfPowerFrequency',2000,'FilterOrder',10,'DesignMethod','butter');
% temp_data = filtfilt(filt,temp_data);
temp_data = downsampleTLab(temp_data,Fs/1000,3); %traditional downsampling --> ~1k
else
    
end
% wheel_raw = unwrapBeh(data);
[wheel,flagmtx ]= unwrapBeh(temp_data(1:end-n_delete));
% [finalData, flagMatrix] = unwrapBeh(rawData);
% a smooth measure of number of periods (usually, wheel rotations) in the signal

% Note:
% theoretically speaking it should be  <  0.4 * newsamplingrate
% Filtering before downsampling is required to avoid digital aliasing

wheel = downsampleTLab(wheel,dsRate/(Fs/1000),2); % binning and mean 'downsample' to desired freq

distance = wheel*circum; %distance on the 
vel = getVel(distance', Fs/(dsRate/(Fs/1000)), winSize); % moving window mean and median filtered velocity trace

    minRest = params.beh.minRestTime * Fs/(dsRate/(Fs/1000));
    minRun = params.beh.minRunTime * Fs/(dsRate/(Fs/1000));
onsets = [];
offsets = [];
% [onsets,offsets] = getOnsetOffset(abs(vel),velThres,minRest,minRun,finalOnset);

return
%% 
figure;
subplot(4,1,1);
plot(data);
subplot(4,1,2);
plot(temp_data);
subplot(4,1,3);
plot(wheel);
subplot(4,1,4);
plot(vel)
%% original from Pratik ... why need to filter the data
 wheel = unwrapBeh(wheel);
    lpFilt = designfilt('lowpassiir','SampleRate',rawFs,'FilterOrder',10,'HalfPowerFrequency',10);
    wheel = filtfilt(lpFilt,wheel);
    wheel = downsample_TLab(wheel,dsRate,dsType);
    wheel = wheel*circum;
    data.final(n).wheel = wheel;
    vel = getVel(wheel,Fs,winSize);
    minRest = params.beh.minRestTime * Fs; minRun = params.beh.minRunTime * Fs;
    [onsets,offsets] = getOnsetOffset(abs(vel),velThres,minRest,minRun,finalOnset);
    data.final(n).vel = vel';
    if ~isfield(data.final(n),'time')
        timeVec = [1:length(vel)];
        data.final(n).time = timeVec'/Fs;
    end        
    data.final(n).beh.onsets = onsets;
    data.final(n).beh.offsets = offsets;