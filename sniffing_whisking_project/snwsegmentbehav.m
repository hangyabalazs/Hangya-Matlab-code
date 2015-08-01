function snwsegmentbehav
%SNWSEGMENTBEHAV   Behavior time stamps.
%   SNWSEGMENTBEHAV adds time stamps from behavioral data to the 'segments'
%   structures created by SNWSNIFFPLOT2. Time of pellet drop and auditory
%   cue onset are added.
%
%   See also SNWSNIFFPLOT2 and SNWSEGMENTSPEED.

% Load all segments
global DATADIR
fn = [DATADIR 'HSW\allsegments_speed.mat'];
load(fn)

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running SNWSEGMENTBEHAV...');  % progress indicator
global WB
WB(end+1) = wb;

% Behavioral data for each segment
NumSeg = length(allsegments);   % number of segments
allsegments(1).pellet = [];  % initialize fields
allsegments(1).tone = [];
for k = 1:NumSeg
    allsegments(k) = segmentbehav(allsegments(k));   % behavioral data for the segment
    waitbar(k/NumSeg)   % update progress indicator
end
close(wb)   % close progress indicator

% Save
fn = [DATADIR 'HSW\allsegments_behav2.mat'];
save(fn,'allsegments')

% -------------------------------------------------------------------------
function cseg = segmentbehav(cseg)

% Data import
global BEHAVDATA
global DATADIR
rrat = cseg.rat;  % rat ID
rdate = cseg.session;   % session date
if ~isempty(BEHAVDATA) && isequal(BEHAVDATA.rat,rrat) && isequal(BEHAVDATA.date,rdate)
    pellet = BEHAVDATA.pellet;  % if the data was loaded in the previous round, use it
    tone = BEHAVDATA.tone;
else
    fn = [DATADIR 'HSW\behav\DATA_' rrat '_' rdate '.mat'];  % if not, load it
    try
        load(fn)
    catch ME
        disp('Failed to load behavioral data.')
        disp(ME.message)
        pellet = NaN;   % load problem: usual cause is missing file
        tone = NaN;
    end
    if exist('valid_pellet_on','var')   % if there is behavioral data
        pellet = valid_pellet_on - time(1);   % time stamp of pellet drop
    else
        pellet = NaN;
    end
    if exist('valid_tone_on','var')   % if there is behavioral data
        tone = valid_tone_on - time(1);   % time stamp of cue onset
    else
        tone = NaN;
    end
end

% Store data as global
BEHAVDATA.pellet = pellet;   % time stamp of pellet drop
BEHAVDATA.tone = tone;   % time stamp of cue onset
BEHAVDATA.rat = rrat;   % animal ID
BEHAVDATA.date = rdate;   % session date

% Add behavioral data to the 'segments' structure
cseg.pellet = pellet;
cseg.tone = tone;