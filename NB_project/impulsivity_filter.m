function [valid_trials, TE] = impulsivity_filter(TE)
%IMPULSIVITY_FILTER   Finds impulsive trials.
%   [VALID_TRIALS, TE] = IMPULSIVITY_FILTER(TE) find trials where at least
%   one ITI was >30s in the past 10 trials ('impulsive' trials). Indices of
%   all other trials (VALID_TRIALS) and trial event structure with
%   impulsive trials removed (TE) are returned.
%
%   See also COND_ACCURACY_FILTERED.

%   Edit log: BH 7/8/11

% Find 'impulsive' trials
NUMtrials = length(TE.TrialStart);
win = 10;
inx = [];
valid_trials = [];
for ibl = 1:NUMtrials-win
    if any(TE.TotalITI(ibl:ibl+win)>30)
        inx = [inx ibl+10];   % find trials where at least one ITI was >30s in the past 10 trials
    else
        valid_trials = [valid_trials ibl+win];
    end
end

% Drop 'impulsive' trials
fnm = fieldnames(TE);
for k = 1:length(fieldnames(TE))
    TE.(fnm{k})(inx) = [];
end