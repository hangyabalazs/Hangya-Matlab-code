function snwtonestruct
%SNWTONESTRUCT   Behavioral data structure.
%   SNWTONESTRUCT restructures behavioral data of sniffing-whisking
%   experiments. It saves a structure containing an entry for each tone
%   presentation and stores related information of animal ID, session data,
%   time stamp, trial validity and trial category based on pre-pellet
%   sniffing frequency.
%
%   See also SNWPELLETSTRUCT, SNWSEGMENT_BEHAV and SNWSEGMENT_SPEED.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection2.xlsx'];
[tbl0 tbl] = xlsread(tblfile);   % read session information

% Load behavior information
global DATADIR
fn = [DATADIR 'HSW\snwdatainfo2a.mat'];
load(fn)

% Call 'pstruct' to make the structure
numSessions = size(tbl,1);   % number of sessions from all rats
tones = struct([]);   % preallocate pellet structure
for iS = 1:numSessions    % loop through sessins
    rrat = tbl{iS,1};   % rat ID
    rdate = tbl{iS,2};  % session date
    sessionstart = tbl0(iS,4);   % session start time stamp
    if isequal(rrat,'P5')   % no behavior for this animal
        continue
    end
    
    sessnames = cellfun(@(s)s(9:end-4),snwdata.sessname,'UniformOutput',false);   % session dates
    sessinx = find(strcmp(rdate,sessnames));   % index for current session session
    tone = snwdata.toneon{sessinx};   % auditory cue tone
    invalid = snwdata.invalid_trials{sessinx};   % invalid trials
    lopresniff = snwdata.lowpresniff{sessinx};   % low sniff frequency before pellet drop
    hipresniff = snwdata.hipresniff{sessinx};   % high sniff frequency before pellet drop
    if ~isempty(tone)   % if pellet drop time stamps are available
        lp = pstruct(rrat,rdate,sessionstart,tone,invalid,lopresniff,hipresniff);
        tones = [tones lp];
    end
end

% Save
fn = [DATADIR 'HSW\tones.mat'];
save(fn,'tones')

% -------------------------------------------------------------------------
function tones = pstruct(rrat,rdate,sessionstart,tone,invalid,lopresniff,hipresniff)

NumTones = length(tone);   % number of tone presentations
tones = struct([]);
for iPD = 1:NumTones    % loop through all pellet drops
    tones(iPD).rat = rrat;   % rat ID
    tones(iPD).session = rdate;  % session date
    tones(iPD).timestamp = tone(iPD) * 1e-6 - sessionstart;   % tone onset time stamps
    tones(iPD).validity = ~ismember(iPD,invalid);   % trial validity
    tones(iPD).lopresniff = ismember(iPD,lopresniff);   % low sniff frequency before pellet drop
    tones(iPD).hipresniff = ismember(iPD,hipresniff);   % high sniff frequency before pellet drop
end