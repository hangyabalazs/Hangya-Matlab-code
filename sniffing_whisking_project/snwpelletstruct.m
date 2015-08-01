function snwpelletstruct
%SNWPELLETSTRUCT   Behavioral data structure.
%   SNWPELLETSTRUCT restructures behavioral data of sniffing-whisking
%   experiments. It saves a structure containing an entry for each pellet
%   drop and stores related information of animal ID, session data, time
%   stamp, trial validity and trial category based on pre-pellet sniffing
%   frequency.
%
%   See also SNWSEGMENT_BEHAV and SNWSEGMENT_SPEED.

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
pelletdrops = struct([]);   % preallocate pellet structure
for iS = 1:numSessions    % loop through sessins
    rrat = tbl{iS,1};   % rat ID
    rdate = tbl{iS,2};  % session date
    sessionstart = tbl0(iS,4);   % session start time stamp
    if isequal(rrat,'P5')   % no behavior for this animal
        continue
    end
    
    sessnames = cellfun(@(s)s(9:end-4),snwdata.sessname,'UniformOutput',false);   % session dates
    sessinx = find(strcmp(rdate,sessnames));   % index for current session session
    pellet = snwdata.pelleton{sessinx};   % pellet drop time stamps
    invalid = snwdata.invalid_trials{sessinx};   % invalid trials
    lopresniff = snwdata.lowpresniff{sessinx};   % low sniff frequency before pellet drop
    hipresniff = snwdata.hipresniff{sessinx};   % high sniff frequency before pellet drop
    if ~isempty(pellet)   % if pellet drop time stamps are available
        lp = pstruct(rrat,rdate,sessionstart,pellet,invalid,lopresniff,hipresniff);
        pelletdrops = [pelletdrops lp];
    end
end

% Save
fn = [DATADIR 'HSW\pelletdrops2.mat'];
save(fn,'pelletdrops')

% -------------------------------------------------------------------------
function pelletdrops = pstruct(rrat,rdate,sessionstart,pellet,invalid,lopresniff,hipresniff)

NumPelletDrops = length(pellet);   % number of pellet drops
pelletdrops = struct([]);
for iPD = 1:NumPelletDrops    % loop through all pellet drops
    pelletdrops(iPD).rat = rrat;   % rat ID
    pelletdrops(iPD).session = rdate;  % session date
    pelletdrops(iPD).timestamp = pellet(iPD) * 1e-6 - sessionstart;   % pellet drop time stamps
    pelletdrops(iPD).validity = ~ismember(iPD,invalid);   % trial validity
    pelletdrops(iPD).lopresniff = ismember(iPD,lopresniff);   % low sniff frequency before pellet drop
    pelletdrops(iPD).hipresniff = ismember(iPD,hipresniff);   % high sniff frequency before pellet drop
end