% function [Info, Result, WS_data] = HL_FP_Process_Lockin_wisos(ses_fn, FP_ch_name, FP_ref_name)
% function to process WS data get dF/F for using locking method imaging
% green and isosbestic signals
% for next step analysis
% it uses several helper function and pipeline Pratik wrote
% now it works for continues imaging session or single sweep imaging session
% Stim_ch_name: for HL rig Blue_ctx or RED_LED
%
%
%   OUTPUT: (need update)
%         Info.
%               
%         Result.ts_ds = ts_ds;
%         Result.df_F_ds = df_F_ds;
%         Result.FP_clean_inpaint = FP_clean_inpaint;
%         Result.df_F = df_F; % without downsampling, F_baseline, FP_filter
%         Result.F_baseline = F_baseline; % F0 without downsampling
%         Result.FP_filter = FP_filter; % raw data after filter
%         Result.Stim_ts = Stim_ts;
%         Result.trial_label = trial_label;
%         Result.trial_type = trial_type;
%         Result.idx_byTrialType = idx_byTrialType;
%   
%         WS_data. 
% 
% 
% Function dependency:
%   HL_FP_loadWS_parseData.m
%   HL_FP_df_cw.m
% status: complete
% Haixin Liu 2019-10
function [Info, Result, WS_data] = HL_FP_Process_Lockin_wisos(ses_fn, FP_ch_name, FP_ref_name)
%% default paramters
if nargin < 3
   FP_ch_name = 'FP1';
%    Stim_ch_name = 'Blue_Ctx';
%    stimulus_length = 500;
% elseif nargin >2 && nargin < 5
%     help HL_FP_Process_Stim_CW
%     error('input number not matched')
   FP_ref_name{1} = 'FP1_ref';
   FP_ref_name{2} = 'FP2_ref';
   
end

fprintf(2,'Channel name used:\nFP ch: %s\n', FP_ch_name);
fprintf(2,'Default Config: FP1_ref drives blue LED, FP2_ref drives UV LED\n');
%% checking plot flags
flag_check_raw = 0;
flag_check_demo = 0;
flag_check_fit_points = 0; 
flag_check_fit = 0;
flag_cmp = 0;
%% load and parse data
ses_data = HL_FP_loadWS_parseData(ses_fn);
% [DATA] = HL_FP_loadWS_parseData(fn) HL_loadWS_parseData_Csti_FP1ch

% get chn idx
idx_FP_ch = find(cellfun(@(x) contains(x,FP_ch_name) && ~contains(x,'ref'), ses_data.ch_names));
idx_FPref_ch = nan(size(FP_ref_name));
for ii = 1:length(FP_ref_name)
idx_FPref_ch(ii) = find(cellfun(@(x) contains(x,FP_ref_name{ii}) , ses_data.ch_names));
end

[WS_trial, ~, map_num_used] = HL_FP_parseWSStiLib(ses_data.StiLib); % HL_WS_parseStiLib
disp('Ch names');disp(ses_data.ch_names);
% idx_Cstim_ch = find(cellfun(@(x) contains(x,Stim_ch_name), ses_data.ch_names));
%     idx_Bitcode_ch = find(cellfun(@(x) contains(x,'Bitcode'), ses_data.ch_names));

% [FP_clean, Thred, n_trial] = HL_FP_CleanStiArtiFact(ses_data.ch_data(:,idx_FP_ch),ses_data.ch_data(:,idx_Cstim_ch),ses_data.sr); % HL_cleanFP_StiArtiF
% interpolate NaNs
% FP_clean_inpaint = inpaint_nans(FP_clean,5); % just use nearest average to fill NaNs

% Use Pratik's pipeline to coordinate across Tlab
% [df_F_ds, ts_ds, df_F, F_baseline, FP_filter] = HL_FP_df_cw (rawFP, ts, rawFs, system_baseline, ...
%                                                                         lpCut, filtOrder, interpType, fitType, winSize, winOv, basePrc)
% [df_F_ds, ts_ds, df_F, F_baseline, FP_filter] = ...
%     HL_FP_df_cw (ses_data.ch_data(:,idx_FP_ch), ses_data.ts, ses_data.sr);
%% check raw singal
if flag_check_raw
    figure; a = [];
    for ii = 1:length(idx_FPref_ch)
    a(ii)= subplot(length(idx_FPref_ch)+1,1,ii);
    plot(ses_data.ts, ses_data.ch_data(:,idx_FPref_ch(ii)));   
    title(['REF FP: ' FP_ref_name{ii}])
    end
    ii = ii +1;
    a(ii)= subplot(length(idx_FPref_ch)+1,1,ii);
    plot(ses_data.ts, ses_data.ch_data(:,idx_FP_ch));
    title('FP raw');
    xlabel('Time (s)')
    linkaxes(a,'x')
end
%%
for ii = 1:length(FP_ref_name)
    if ii == 1
        [df_F_ds{ii}, ts_ds, df_F{ii}, F_baseline{ii}, rawFP_demod_filtered{ii}, params] = ...
            HL_FP_df_lockin (ses_data.ch_data(:,idx_FP_ch), ses_data.ch_data(:,idx_FPref_ch(ii)),  ses_data.ts, ses_data.sr);
    else
        [df_F_ds{ii}, ~, df_F{ii}, F_baseline{ii}, rawFP_demod_filtered{ii}] = ...
            HL_FP_df_lockin (ses_data.ch_data(:,idx_FP_ch), ses_data.ch_data(:,idx_FPref_ch(ii)),  ses_data.ts, ses_data.sr);
    end
end
%% check demodulated filtered FP signal
if flag_check_demo
    figure; a = [];
    for ii = 1:length(df_F_ds)
    a(ii)= subplot(length(df_F_ds)+1,1,ii);
    plot(ses_data.ts, rawFP_demod_filtered{ii});   
    title(['demod&filtered FP : ' FP_ref_name{ii}])
    end
    
    xlabel('Time (s)')
    linkaxes(a,'x')
    
    a(1+length(df_F_ds))= subplot(length(df_F_ds)+1,1,1+length(df_F_ds));
    hold on;
    for ii = 1:length(df_F_ds)
    
    plot(ts_ds, df_F_ds{ii});   
%     title(['demod&filtered FP dF/F: ' FP_ref_name{ii}])
    end
%     legend(num2str(1:length(df_F_ds)))
    title('dF/F calc using sliding window method');
    ylabel('dF/F (%)')
    xlabel('Time (s)')
    linkaxes(a,'x')
end
%% iso method => get baseline using fitting
% ? negatively correlated ? 
idx_fit = ~isnan(rawFP_demod_filtered{2} + rawFP_demod_filtered{1});

if flag_check_fit_points
figure;
histogram2(rawFP_demod_filtered{2}, rawFP_demod_filtered{1},'DisplayStyle','tile')
pbaspect([1 1 1]); hold on;

[fitobject,gof] = fit(rawFP_demod_filtered{2}(idx_fit),rawFP_demod_filtered{1}(idx_fit),'poly1');
plot([min(rawFP_demod_filtered{2}) max(rawFP_demod_filtered{2})],...
    fitobject([min(rawFP_demod_filtered{2}) max(rawFP_demod_filtered{2})]), '-b');

[fitobject,gof] = fit(rawFP_demod_filtered{2}(idx_fit),rawFP_demod_filtered{1}(idx_fit),'poly1','Robust','Bisquare');
plot([min(rawFP_demod_filtered{2}) max(rawFP_demod_filtered{2})],...
    fitobject([min(rawFP_demod_filtered{2}) max(rawFP_demod_filtered{2})]), '-r');
title('Robust:b-OFF, r-Bisquare')
xlabel('FP Ref 2: iso (UV)');
ylabel('FP Ref 1: green');

end

[fitobject,gof] = fit(rawFP_demod_filtered{2}(idx_fit),rawFP_demod_filtered{1}(idx_fit),'poly1');

rawFP_demod_filtered_fit = fitobject(rawFP_demod_filtered{2});
%% dF/F calc 
df_F_demod = ((rawFP_demod_filtered{1} - rawFP_demod_filtered_fit) ./ rawFP_demod_filtered_fit)*100;
% downsample
dsRate = ses_data.sr/params.FP.ds2;
df_F_demod_ds =  downsample(df_F_demod,dsRate);

%%
if flag_check_fit
    figure; 
    a= [];
    a(1)=subplot(3,1,1);
    plot(ses_data.ts, rawFP_demod_filtered{1}, 'k');  
    hold on;
    plot(ses_data.ts, rawFP_demod_filtered_fit, 'b');
    legend('FP1','FP2 fit');
    title('Demodu & filtered')
    a(2) = subplot(3,1,2);
    plot(ses_data.ts, rawFP_demod_filtered{2}, 'k');
     hold on;
    plot(ses_data.ts, rawFP_demod_filtered_fit, 'b');
    legend('FP2', 'FP2 fit');
    a(3) = subplot(3,1,3);
    plot(ses_data.ts, df_F_demod, 'g');
    legend('dF/F');
    xlabel('Time (s)');
    ylabel('dF/F (%)')
    
    linkaxes(a,'x');
   
end
%% 
if flag_cmp
    figure; a= [];
    a(1)=subplot(2,1,1);
    hold on;    
    plot(ses_data.ts, df_F{1}, 'b');
    plot(ses_data.ts, df_F_demod, 'r');
    
    legend('moving window','iso fit');
%     xlabel('Time (s)');
    ylabel('dF/F (%)')
    
    a(2)=subplot(2,1,2);
    plot(ses_data.ts, df_F{1} - df_F_demod, 'k');
        xlabel('Time (s)');
    ylabel('Diff in dF/F (%)');
    title('mov.win. - fit');
        linkaxes(a,'x');

end
%%

% get trial onset of BitCode ...
% TrialNumber = HL_ReadBitCode(ses_data.ch_data(:,idx_Bitcode_ch),ses_data.sr);
% [Stim_ts, ~]= HL_FP_ParseSti(ses_data.ch_data(:,idx_Cstim_ch),ses_data.sr, Thred, stimulus_length);

% [trial_label, trial_type, idx_byTrialType] = HL_FP_GenerateTrialTypeIndex(size(Stim_ts,1), WS_trial);

% record data


%% return useful result
Info.fitobject = fitobject;

% Info.auto_baseline =auto_baseline ;
% Info.sys_noise_est =sys_noise_est ;
% Info.WS_trial = WS_trial;
% Info.Stim_Thred = Thred;
% Info.FP_trial_window_length = FP_trial_window_length;
% Info.FP_trial_window_start = FP_trial_window_start;
% Info. = ;
% Info. = ;
% Info. = ;
% Info. = ;

Result.df_F_demod_ds = df_F_demod_ds;
Result.df_F_demod = df_F_demod;
Result.ts_ds = ts_ds;
Result.ts = ses_data.ts;
Result.rawFP_demod_filtered = rawFP_demod_filtered; % raw data after filter
Result.rawFP_demod_filtered_fit = rawFP_demod_filtered_fit;
%
Result.df_F_ds = df_F_ds;
% Result.FP_clean_inpaint = ses_data.ch_data(:,idx_FP_ch); % just use raw to return % FP_clean_inpaint;
Result.df_F = df_F; % without downsampling, F_baseline, FP_filter
Result.F_baseline = F_baseline; % F0 without downsampling
%
Result.params = params;

% Result.Stim_ts = Stim_ts;
% Result.trial_label = trial_label;
% Result.trial_type = trial_type;
% Result.idx_byTrialType = idx_byTrialType;
% Result.n_trial = n_trial;
% Result.Trial_FP_15s = Trial_FP_15s;
% Result.Trial_FP_1m = Trial_FP_1m;
% Result.FP_x_plot = FP_x_plot;


WS_data = ses_data;
