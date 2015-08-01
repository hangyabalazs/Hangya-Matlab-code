function valid_trials = filterTrials(cellid,varargin)
%FILTERTRIALS   Filtering trials.
%   VALID_TRIALS = FILTERTRIALS(CELLID,PARAM,VALUE) filters the events for
%   a given cell (CELLID) based on the following parameter, value pairs
%   (default values are given below):
%       'event_type', 'stim' - select whether which event type ('stim' or
%           'trial') will be used
%       'event', 'PulseOn' - reference event
%       'event_filter', 'none' - filter name; see implemented filter types 
%           below; intersection of filters is used for a cell array of
%           filter names
%       'filterinput', [] - some filters require additional input - see
%           below; filter input has to be given as a struct, with
%           appropriate fields
%
%   Implemented event filters for 'stim' events:
%       'BurstNPulse' - only light pulse trains with a specified number of
%           pulses are used; this number has to be specified as a
%           'BurstNPulse', N (numeric) parameter, value pair given as a
%           struct (s.BurstNPulse = 1)
%       'BurstNPulse_maxPower' - the same as 'BurstNPulse' but light bursts
%           are further restricted to pulses with maximal power applied
%       'minNPulse_maxPower' - the light bursts with minimal number of
%           pulses and maximal power are used
%
%   Implemented event filters for 'trial' events:
%       'animalactive' - filters for minimal percentage of hits (40%=4 hits) 
%           over a 10 go-trial moving window
%       'fastGo' - upper percentile of GoRT; 'filterinput' specifies the
%           percentile threshold
%       'slowGo' - lower percentile of GoRT; 'filterinput' specifies the
%           percentile threshold
%
%   'custom' filter:
%       Set event_filter to 'custom' and use the 'filterinput' argument to
%       pass a logical expression for filtering trials (e.g.
%       TE.Hit==1&TE.GoRT<0.5) - see SELECTTRIAL for details.
%
%   See also PARSEVALIDTRIALS and SELECTTRIAL.

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'event','PulseOn',@ischar)   % reference event
addParamValue(prs,'event_type','stim',@ischar)   % event type ('stim' or 'trial')
addParamValue(prs,'event_filter','none',@(s)ischar(s)|iscellstr(s))   % event filter
addParamValue(prs,'filterinput',[])   % some filters need additional input
parse(prs,cellid,varargin{:})
g = prs.Results;

% Load event structure
event_type = lower(g.event_type(1:4));
switch event_type
    case 'stim'
        
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');
        catch ME
            disp('There was no stim protocol for the session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');
        catch ME
            disp('There was no behavioral protocol for ths session.')
            error(ME.message)
        end
        
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

% Filter events
if iscellstr(g.event_filter)
    fld = fieldnames(VE);
    valid_trials = 1:length(VE.(fld{1}));
    for flts = 1:length(g.event_filter)
        valid_trials = intersect(valid_trials,...
            main(VE,g.event,g.event_filter{flts},g.filterinput{flts}));
    end
else
    valid_trials = main(VE,g.event,g.event_filter,g.filterinput);
end

% -------------------------------------------------------------------------
function valid_trials = main(VE,event,event_filter,filterinput)

% Filter trials
switch event_filter
    case 'minNPulse_maxPower'
        minfreq = min([VE.BurstNPulse]);
        maxpow = max([VE.PulsePower]);
        inx = ~isnan(VE.(event)) & VE.BurstNPulse==minfreq & VE.PulsePower==maxpow;
        valid_trials = find(inx);
    case 'BurstNPulse_maxPower'
        maxpow = max([VE.PulsePower]);
        inx = ~isnan(VE.(event)) & VE.BurstNPulse==filterinput.BurstNPulse & VE.PulsePower==maxpow;
        valid_trials = find(inx);
    case 'BurstNPulse'
        inx = ~isnan(VE.(event)) & VE.BurstNPulse==filterinput.BurstNPulse;
        valid_trials = find(inx);
    case 'animalactive'
        % filter for TrialEvents; filters for minimal percentage of hits
        % (40%=4 hits) over a 10 go-trial moving window
        wn = 10;  % integration window
        hitinx = find(VE.Hit==1);  % indices for Hits
        missinx = find(VE.Miss==1);   % indices for False Alarms
        hminx = sort([hitinx missinx]);   % indices for Hits or False Alarms
        NumTrials = length(VE.Hit);
        pchit = nan(1,NumTrials);
        for t = 1:NumTrials
            ahinx = hminx(hminx<=t);  % actual window
            if length(ahinx) >= wn
                cahinx = ahinx(end-wn+1:end);
            else   % in the first wn hits, use (partially) forward window
                cahinx = hminx([find(hminx<=t) find(hminx>t,wn-length(ahinx),'first')]);
            end
            pchit(t) = length(intersect(cahinx,hitinx)) / wn;  % percent Hit / go tones
        end
        valid_trials = find(pchit>0.4);
    case 'fastGo'
        prcthres = filterinput;   % percentile threshold
        rts = VE.GoRT;
        rts = rts(~isnan(rts));   % all RT values
        rtthres = prctile(rts,prcthres*100);   % RT threshold
        valid_trials = find(VE.GoRT>=rtthres(1)&VE.GoRT<=rhthres(2));   % RTs of the percentile interval
    case 'slowGo'
        prcthres = filterinput;   % percentile threshold
        rts = VE.GoRT;
        rts = rts(~isnan(rts));   % all RT values
        rtthres = prctile(rts,prcthres);   % RT threshold
        valid_trials = find(VE.GoRT<rtthres);   % RTs of the lower percentile        
    case 'none'
        valid_trials = find(~isnan(VE.(event)));
    case 'custom'
        valid_trials = selecttrial(VE,filterinput);
    otherwise
        error('Unknown event filter.')
end