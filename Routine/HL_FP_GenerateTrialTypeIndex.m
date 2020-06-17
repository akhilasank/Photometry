% function [trial_label, trial_type, idx] = HL_FP_GenerateTrialTypeIndex(n_trial, trial)
% generate the trail type numbers
% INPUT: 
%  n_trial: total trial number of this session;
%  trial: -struct, output from HL_WS_parseStiLib.m
%             trial.type: cell arrays of used stimulus map
%             tiral.label: idx number in .type of a full stimulus sequence
%
% Haixin Liu 2019-09

function [trial_label, trial_type, idx] = HL_FP_GenerateTrialTypeIndex(n_trial, trial)

%%%% If number of inputs is less than 2 (only input was number of trials), all trials have the same trial type and are labelled as 1.
if nargin < 2
    % just n_trial, make the same type for all
    trial_label = 1;
    trial_type{1} = 'One';
    idx.(trial_type{1}) = 1:n_trial;
else
% fix names to put into a structure
%%%% Otherwise there are multiple trial types.
%%%% Change the '-' in trial type name to '_'.
for ii = 1:length(trial.type)
    if any(strfind(trial.type{ii}, '-'))
        trial.type{ii}(strfind(trial.type{ii}, '-')) = '_';
    end
end
%%%%%%%%%%%%%
%%%% Using the length of the trial label, identify the number of block by dividing th number of total trials by the number of labels. Create an 1xn_block array of repeated copies of 'trial.label' (the total number of elements should equal the number of trials). @HaixinLiuNeuro is this correct?
trail_label = trial.label;
trial_type = trial.type;
n_block = ceil(n_trial/length(trail_label));
trial_label = repmat(trail_label,1,n_block);
trial_label = trial_label(1:n_trial);

% clean trial types, as some cannot be structure name
%%%% change the '.' in trial type name to '_'.
%%%% change the '@' in trial type name to 'a'.
%%%% Index the trial type to the trial label. 
for ii = 1:length(trial_type)
    if any(strfind(trial_type{ii}, '.'))
        trial_type{ii}(strfind(trial_type{ii}, '.')) = '_';
    end
    if any(strfind(trial_type{ii}, '@'))
        trial_type{ii}(strfind(trial_type{ii}, '@')) = 'a';
    end
    idx.(trial_type{ii}) = find(trial_label==ii);
end
end    
