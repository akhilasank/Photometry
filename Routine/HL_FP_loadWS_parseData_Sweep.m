% function [DATA] = HL_FP_loadWS_parseData_Sweep(fn)
% load wavesurfer data using WS matlab package (ws.loadDataFile)
% OUTPUT: DATA: -struct
%           ch_names: N�1 cell for recording channels
%           sr: sampling rate, sample per second
%           StiLib: [1�1 struct] copy of stimulus library for figureing out
%                   stimulus protocol later
%           ch_data: -matrix data point X channel 
%           ts: time stamp in second, a vector data point X 1. time
%               starting from 0. 
% INPUT: fn: full file name of the .h file 
%
% Haixin Liu 2019-12
% revised from HL_FP_loadWS_parseData.m (for continuous recording which is just one sweep)
% to read in data with the recording config as sweep
% .ts start from 0, since every sweep is the same just use one
% also return the clock time for each sweep start to use later if needed

%% 
function [DATA] = HL_FP_loadWS_parseData_Sweep(fn)
% warning('This function currently only deal with single sweep/continuous recording')
% Use WS package to read in
disp('Reading WS file ...')
tic
s = ws.loadDataFile(fn);
toc
% Process for specific needs

DATA.ch_names = s.header.AIChannelNames(find(s.header.IsAIChannelActive));

% s.header.AOChannelNames

DATA.sr = s.header.AcquisitionSampleRate;

if s.header.AreSweepsContinuous
    disp('Continuous recording');
else
    disp('Sweep-based recording')
end
% s.header.ClockAtRunStart : [yyyy mm dd hh mm ss.ffff]
% s.header.StimulationSampleRate

%stimulus used, need a easy way to parse out, but now can just grab it 
DATA.StiLib = s.header.StimulusLibrary;

% check sweep #
temp_fieldn = fieldnames(s);
temp_idx = find(cellfun(@(x) ~isempty(strfind(x, 'sweep')), temp_fieldn));

% if length(temp_idx) > 1 % for now dealing with single sweep session: continues recording
%     help HL_FP_loadWS_parseData
%     warning('multiple sweeps!!! ONLY return the FIRST sweep!!!')
%     temp_idx = temp_idx(1);
% end
for ii = 1:length(temp_idx)
DATA.ch_data{ii} = s.(temp_fieldn{temp_idx(ii)}).analogScans; % data points x channel 
DATA.sweep_start_ts(ii) = s.(temp_fieldn{temp_idx(ii)}).timestamp; % in sec from recording start
end
%define time start from time 0? first time point is 0, 2nd is
%0+1/sample_rate
DATA.ts = [0:1:(size(s.(temp_fieldn{temp_idx(1)}).analogScans,1)-1)]'/s.header.AcquisitionSampleRate;
disp('Define time start at 0 s')

%% HL 2019-9-4 add in module to read digital input channel as well
% warning('Digigal INPUT not included yet'); return
%{
 from Adam: 
It multiplexes the two inputs into a single uint8 per sample.
 The lowest-order but should be the first signal, 
the next-lowest should be the second signal, 
if I remember correctly.  
In matlab, you can use the bitget() function to extract individual bits.
%}
if ~isempty(s.header.DIChannelNames)
    disp('Having DI channels, read in [only support 8 DI channels for now]')
    warning('Not implemented yet, RETURNED'); return
% s.sweep_0001.digitalScans;
% s.header.DIChannelNames;
tic
DI_convert = arrayfun(@(x) bitget(x,8:-1:1,'uint8')', s.(temp_fieldn{temp_idx}).digitalScans, 'UniformOutput',false);
% bitget(s.sweep_0001.digitalScans(1),8:-1:1,'uint8')
DI_convert_temp = cat(2,DI_convert{:});
toc
% append to ch data and names
for i_DI = 1:length(s.header.DIChannelNames)
    DATA.ch_data = cat(2,DATA.ch_data, double(DI_convert_temp(end-i_DI+1,:)')); % IMPORTANT change to double format
    DATA.ch_names = cat(1, DATA.ch_names, s.header.DIChannelNames(i_DI));
end
else
    
end
%%


return
%% old stuff
%{
figure;
a(1) = subplot(4,1,1); 
plot(s.(temp_fieldn{temp_idx}).analogScans(1:50*2000,1));
title('Camera Frame')
a(2) = subplot(4,1,2); 
plot(s.(temp_fieldn{temp_idx}).analogScans(1:50*2000,2));
title('Ctx Sti')
a(3) = subplot(4,1,3); 
plot(s.(temp_fieldn{temp_idx}).analogScans(1:50*2000,3));
title('Bpod BitCode')
a(4) = subplot(4,1,4); 
plot(s.(temp_fieldn{temp_idx}).analogScans(1:50*2000,4));
title('Photometry')
ylabel('Volt')
xlabel('Data points')
linkaxes(a,'x');
%% 
figure;
a(1) = subplot(2,1,1); 
plot(DATA.ts, FP_clean);
hold on; 
plot(DATA.ts, to_exclude*4);
title('Cleaned FP')
a(2) = subplot(2,1,2); 
plot(DATA.ts, Sti);
title('Ctx Sti')
linkaxes(a,'x');




%}
