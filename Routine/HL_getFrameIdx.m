function [frame_idx] = HL_getFrameIdx(Frametimes, event_ts)
% function [frame_idx] = getFrameIdx(Frametimes, event_ts)
%
% use time stamp and frame time stamps to get the nearest frame for each
% event
% return idx in Input Frametimes
%%
for ii = 1:size(event_ts,1)
    [~,idx_temp]=min(abs(Frametimes - event_ts(ii,1)));
    frame_idx(ii,:) = idx_temp(1);
    
end