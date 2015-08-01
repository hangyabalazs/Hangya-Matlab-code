function firingrate_call
%FIRINGRATE_CALL   Call FIRINGRATE.
%   FIRINGRATE_CALL calls FIRINGRATE for cholinergic neurons.
%
%   See also FIRINGRATE.

% Cholinergic neurons
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells (n = 22)
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes

% Call FIRINGRATE
NumChAT = length(ChAT);   % number of cells
[hit_fr fa_fr miss_fr cr_fr] = deal(nan(NumChAT,2));
prewindow = [-0.5 0];   % windows
postwindow = [0 0.5];
align_event = 'StimulusOn';   % reference event
for iC = 1:NumChAT
    cellid = ChAT{iC};
    disp(cellid)
%     align_event = findAlignEvent_negfeedback_gonogo(cellid);
    fr = firingrate(cellid,'trial',align_event,prewindow,postwindow,...
        'event_filter','custom','filterinput','Hit==1');
    hit_fr(iC,1:2) = nanmean(fr);
    fr = firingrate(cellid,'trial',align_event,prewindow,postwindow,...
        'event_filter','custom','filterinput','FalseAlarm==1');
    fa_fr(iC,1:2) = nanmean(fr);
    fr = firingrate(cellid,'trial',align_event,prewindow,postwindow,...
        'event_filter','custom','filterinput','Miss==1');
    miss_fr(iC,1:2) = nanmean(fr);
    fr = firingrate(cellid,'trial',align_event,prewindow,postwindow,...
        'event_filter','custom','filterinput','CorrectRejection==1');
    cr_fr(iC,1:2) = nanmean(fr);
end

% Statistics
keyboard
boxstat(hit_fr(:,1),hit_fr(:,2),'pre','post')
boxstat(fa_fr(:,1),fa_fr(:,2),'pre','post')
boxstat(cr_fr(:,1),cr_fr(:,2),'pre','post')
boxstat(miss_fr(:,1),miss_fr(:,2),'pre','post')