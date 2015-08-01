function snwsegmentspeed
%SNWSEGMENTSPEED   Animal speed.
%   SNWSEGMENTSPEED estimates and saves speed vector for all 'segments'
%   structures created by SNWSNIFFPLOT2. X and Y position data (sampled at
%   1 kHz) is smoothed by a 51 ms moving average and used for speed vector
%   calculation. The speed vectors are added to the 'segments' structures
%   as 'speed' fields.
%
%   See also SNWSNIFFPLOT2 and SNWCONVERTDATA.

% Load all segments
global DATADIR
fn = [DATADIR 'HSW\allsegments.mat'];
load(fn)
sr = 1000;   % sampling rate

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running SNWSEGMENTSPEED...');  % progress indicator
global WB
WB(end+1) = wb;

% Speed for each segment
NumSeg = length(allsegments);   % number of segments
for k = 1:NumSeg
    if ~ismember(allsegments(k).rat,{'R1','P5'})   % no speed data for these rats
        allsegments(k) = segmentspeed(allsegments(k),sr);   % speed for the segment
    else
        allsegments(k).speed = NaN;   % pad with NaN if no speed data is available
    end
    waitbar(k/NumSeg)   % update progress indicator
end
close(wb)   % close progress indicator

% Save
fn = [DATADIR 'HSW\allsegments_speed.mat'];
% save(fn,'allsegments')

% -------------------------------------------------------------------------
function cseg = segmentspeed(cseg,sr)

% Data import
global POSITIONDATA
global DATADIR
rrat = cseg.rat;  % rat ID
rdate = cseg.session;   % session date
if ~isempty(POSITIONDATA) && isequal(POSITIONDATA.rat,rrat) && isequal(POSITIONDATA.date,rdate)
    xpos = POSITIONDATA.xpos;  % if the data was loaded in the previous round, use it
    ypos = POSITIONDATA.ypos;
else
    fn = [DATADIR 'HSW\mat2\DATA_' rrat '_' rdate '.mat'];  % if not, load it
    load(fn)
end

% Store data as global
POSITIONDATA.xpos = xpos;   % x position
POSITIONDATA.ypos = ypos;   % y poistion
POSITIONDATA.rat = rrat;   % animal ID
POSITIONDATA.date = rdate;   % session date

% Check for validity of position data
if ~any(xpos.*ypos)
    cseg.speed = NaN;   % return NaNs if no position data
    return
end

% Check if position data covers the range of the segment
if cseg.end > length(xpos) / sr
    cseg.speed = NaN;   % return NaNs if no position data
    return
end

% Data window
mrg = 1;   % margin
wn = [cseg.start-mrg cseg.end+mrg];  % from segment start to end, add 1s margin
segstinx = max(round(wn(1)*sr),1);  % index for segment start
segndinx = min(round(wn(2)*sr),length(xpos));   % index for segment end
mrg = [cseg.start-segstinx/sr segndinx/sr-cseg.end];   % actual margin

% Speed
mdp = segstinx + (segndinx - segstinx) / 2;  % mid point
cxp = xpos(segstinx:segndinx);  % x position vector for the segment
cyp = ypos(segstinx:segndinx);   % y position vector for the segment
scxp = smooth(cxp,'linear',51);   % smoothed x positions
scyp = smooth(cyp,'linear',51);   % smoothed y positions
dx = diff(scxp);   % derivative of x positions
dy = diff(scyp);   % derivative of y positions
spd = sqrt(dx.^2+dy.^2) * sr;   % speed

% Remove margin
spd = spd(round(mrg(1)*sr)+1:end-round(mrg(2)*sr));
cseg.speed = spd;