function vipbehavcorr_fr
%VIPBEHAVCORR_FR   Behavioral responses of VIP neurons.
%   VIPBEHAVCORR_FR calculates and saves firing rates between the cue tone
%   and the response (feedback) and after the feedback (same window size).
%
%   See also VIPRESPONSESORTER and VIPBEHAVCORR_SCATTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-May-2013

% Result directory
global DATAPATH
resdir = [DATAPATH 'VIP\FRs\'];

% First CellBase: air-puff
choosecb('VIP_gonogo')
loadcb
cellids = CELLIDLIST;

% Align to cue tone
event1 = 'StimOnset';
event2 = 'Lick';

% Firing rates
eventfilter = 'ResponseType==1';  % for Hit, VIP_gonogo CellBase
[fr_Hit_tone1 fr_Hit_fb1] = main(cellids,event1,event2,eventfilter);
eventfilter = 'ResponseType==2';  % for FA, VIP_gonogo CellBase
[fr_FA_tone1 fr_FA_fb1] = main(cellids,event1,event2,eventfilter);

% Second CellBase: foot-shock
choosecb('VIP_gonogo2')
loadcb
cellids = CELLIDLIST;

% Align to cue tone
event1 = 'StimulusOn';
event2 = {'LeftWaterValveOn' 'LeftPortIn'};

% Firing rates
eventfilter = 'ResponseType==1';  % for Hit, VIP_gonogo2 CellBase
[fr_Hit_tone2 fr_Hit_fb2] = main(cellids,event1,event2{1},eventfilter);
eventfilter = 'ResponseType==2';  % for FA, VIP_gonogo2 CellBase
[fr_FA_tone2 fr_FA_fb2] = main(cellids,event1,event2{2},eventfilter);

% Save
FR_Hit_tone = [fr_Hit_tone1 fr_Hit_tone2];  %#ok<NASGU>  % FR between tone and feedback, Hit trials
FR_Hit_fb = [fr_Hit_fb1 fr_Hit_fb2];  %#ok<NASGU>  % FR after feedback, Hit trials
FR_FA_tone = [fr_FA_tone1 fr_FA_tone2]; %#ok<NASGU>  % FR between tone and feedback, FA trials
FR_FA_fb = [fr_FA_fb1 fr_FA_fb2];  %#ok<NASGU>  % FR after feedback, FA trials
fnm = [resdir 'tone_and_fb_rates3.mat'];
save(fnm,'FR_Hit_tone','FR_FA_tone',...
    'fr_Hit_tone1', 'fr_Hit_tone2',...
    'fr_Hit_fb1', 'fr_Hit_fb2',...
    'fr_FA_tone1', 'fr_FA_tone2',...
    'fr_FA_fb1', 'fr_FA_fb2')

% -------------------------------------------------------------------------
function [fr1 fr2] = main(cellids,event1,event2,eventfilter)

% Firing rates
NumCells = length(cellids);
[fr1 fr2] = deal(nan(1,NumCells));
for iC = 1:NumCells
    cellid = cellids{iC};   % cell ID
    disp(cellid)
    
    % Load prealigned spikes
    TE = loadcb(cellid,'TrialEvents');   % load events
    TS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
    event_pos1 = findcellstr(TS.events(:,1),event1);
    event_pos2 = findcellstr(TS.events(:,1),event2);
    stimes1 = TS.event_stimes{event_pos1};   % spike times relative to event1
    stimes2 = TS.event_stimes{event_pos2};   % spike times relative to event2
    event2_times = trialevents2relativetime(TE,event1,event2);   % second event times relative to first event
    valid_trials = filterTrials(cellid,'event_type','trial','event',event1,...
        'event_filter','custom','filterinput',eventfilter);
    NumTrials = length(valid_trials);
    
    % Firing rates
    [spno1 spno2 len lenf1 lenf2] = deal(nan(1,NumTrials));
    for iT = 1:NumTrials
        lst1 = stimes1{valid_trials(iT)};  % spike times in the trial rel. to event1
        lst2 = stimes2{valid_trials(iT)};  % spike times in the trial rel. to event2
        len(iT) = event2_times(valid_trials(iT));  % time between event1 and event2 in the trial
        lenf2(iT) = 0.15;
        lenf1(iT) = 0.05;
        spno1(iT) = length(lst1(lst1>0&lst1<len(iT)));  % number of spikes after event1 within 'len'
        spno2(iT) = length(lst2(lst2>lenf1(iT)&lst2<lenf2(iT)));  % number of spikes after event2 within 'len'
    end
    fr1(iC) = sum(spno1) / sum(len);   % firing rate after first event
    fr2(iC) = sum(spno2) / sum(lenf2-lenf1);   % firing rate after second event 
end