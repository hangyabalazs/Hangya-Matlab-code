function psych_call_newdata
%PSYCH_CALL   Call PSYCH_GONOGO.
%   PSYCH_CALL calls PSYCH_GONOGO to calculate average performance and
%   reaction time in the auditory go/no-go task for all session containing
%   NB cells. See PSYCH_GONOGO for details.
%
%   See also PSYCH_GONOGO.

% Load cellbase
global CELLIDLIST ANALYSES TheMatrix %#ok<NUSED>
if isempty(CELLIDLIST)
    load(getpref('cellbase','fname'));
end

% NB neurons
Lratio = getvalue('Lr_PC');
ID = getvalue('ID_PC');
vldty = getvalue('validity');
area1 = getvalue('Area1');
area2 = getvalue('Area2');
session_type = getvalue('session_type');
inb = isnb(area1,area2);   % select NB cells
% inb = ishdb(area1,area2);   % select NB cells
ptinx = vldty == 1 & ID > 20 & Lratio < 0.15 & inb & ...
    (strcmp(session_type,'feedbackdelay') | strcmp(session_type,'gonogo'));   % good clusters from NB
I = ptinx;
cellids = CELLIDLIST(I);

% Performance
% psych_gonogo2(cellids)
% psych_gonogo_learningcurve(cellids)
psych_gonogo_perftest(cellids)
% psych_gonogo_updating(cellids)
% psych_gonogo_updating_noresp(cellids)
% psych_gonogo_updating_noresp2(cellids)
% psych_gonogo_restarts(cellids)

% -------------------------------------------------------------------------
function I = isnb(area1,area2)

nbareas = {'GP','GP/SI','SI','IC','RT/IC','EP','EA','EAC'};   % areas considered part of the basal forebrain (BF)
I = ismember(area1,nbareas);  % if 'primary' area is part of the BF

% -------------------------------------------------------------------------
function I = ishdb(area1,area2)

nbareas = {'HDB','SI','VP'};   % areas considered part of the basal forebrain (BF)
I = ismember(area1,nbareas);  % if 'primary' area is part of the BF