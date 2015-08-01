function snwmodebehav
%SNWMODEBEHAV   Coupling mode statistics.
%   SNWMODEBEHAV quantifies the clustering of sniffing-whisking coupling
%   modes around behavioral events (pellet drops). It evaluates statistical
%   significance from the point of view of the modes, i.e. is a certain
%   mode more prevalent in a behavioral window normalized to its overall
%   presence than other modes. The overlap between segments of a given mode
%   (see SNWSNIFFPLOT2) with windows relative to behavioral time stamps is
%   assessed. The number of segments overlapping with the behav. windows
%   relative to all segments of that mode is computed. Chi-square test is
%   used for statistical comparison.
%
%   See also SNWSNIFFPLOT2, SNWMODESWITCH and SNWMODESWITCH2.

% Sampling rate
sr = 1000;

% Load all segments
global DATADIR
fn = [DATADIR 'HSW\allsegments_behav.mat'];
load(fn)
fn = [DATADIR 'HSW\switches.mat'];
load(fn)

% Mode indices
gr = [allsegments.group];   % group codes corresponding to all segments
gr05 = find(gr==0.5);   % indices for 1:2 mode
gr1 = find(gr==1);   % indices for 1:1 mode
gr2 = find(gr==2);   % indices for 2:1 mode

% Mid points
stp = [allsegments.start];   % start points
enp = [allsegments.end];   % end points
mip = (stp + enp) / 2;   % middle points

% Behavioral time stamps
pellet = {allsegments.pellet};   % pellet drop
tone = {allsegments.tone};   % auditory cue

% Quantify grouping of modes around pellet drop 
NumSeg = length(allsegments);   % aggregate number of segments
spc = nan(1,NumSeg);   % preallocate
for iS = 1:NumSeg   % loop through all segments
    pelletwin = [pellet{iS}+2; pellet{iS}+8];   % window around all pellet drops
    segwin = [stp(iS); enp(iS)];   % first row: start points; second row: end points
    if ~all(isnan(pelletwin(:)))
        spc(iS) = winoverlap(segwin,pelletwin);   % percent of the segment in the pellet window
    end
end

% Statistics
spc05 = spc(gr05);  % percent overlap for 1:2 mode
spc05(isnan(spc05)) = [];
spc1 = spc(gr1);   % percent overlap for 1:1 mode
spc1(isnan(spc1)) = [];
spc2 = spc(gr2);   % percent overlap for 2:1 mode
spc2(isnan(spc2)) = [];

ranksum(spc05,spc1)   % Mann-Whitney-Wilcoxon tests
ranksum(spc1,spc2)
ranksum(spc05,spc2)

[h p] = b_chi2test2([sum(spc2>0) sum(spc2<=0)],[sum(spc05>0) sum(spc05<=0)],0.05)   % chi square test
[h p] = b_chi2test2([sum(spc1>0) sum(spc1<=0)],[sum(spc05>0) sum(spc05<=0)],0.05)
[h p] = b_chi2test2([sum(spc1>0) sum(spc1<=0)],[sum(spc2>0) sum(spc2<=0)],0.05)

keyboard

% Find close behavioral events
NumSeg = length(allsegments);   % aggregate number of segments
[mpellet1 mpellet2 pellet_dist1 pellet_dist2] = deal(nan(1,NumSeg));   % preallocate
for iS = 1:NumSeg   % loop through all segments
%     [ct mpinx] = min(abs(pellet{iS}-mip(iS)));   % find the closest pellet drop
    mpinx = find(mip(iS)-pellet{iS}>0,1,'last');   % find the closest pellet drop
    if ~isempty(mpinx)   % if there behavioral time stamps are available
        mpellet1(iS) = pellet{iS}(mpinx);   % time stamp of the last pellet drop
        ct = mip(iS) - mpellet1(iS);
        pellet_dist1(iS) = ct;   % time elapsed until the mid point of the segment since the last pellet drop
        
        if mpinx <= length(pellet{iS}) - 1   % if there was an event after the segment
            mpellet2(iS) = pellet{iS}(mpinx+1);   % time stamp of the next pellet drop after the segment
            ct = mip(iS) - mpellet2(iS);
            pellet_dist2(iS) = ct;   % time elapsed until the next pellet drop since the mid point of the segment
        end
    end
end

% -------------------------------------------------------------------------
function [spc ppc] = winoverlap(segwin,pelletwin)

% Calculate length covered by the segments
allwin = [pelletwin segwin];   % "union" of segments
lpw = lebesgue(pelletwin);  % pellet windows
lsw = lebesgue(segwin);  % current segment
law = lebesgue(allwin);   % all segments ("union")

% Calculate overlap
olp = lpw + lsw - law;   % calculate intersect from union

%  Relative overlap
spc = olp / lsw;   % overlap relative to the current segment
ppc = olp / lpw;   % overlap relative to the pellet windows

% -------------------------------------------------------------------------
function len = lebesgue(win)

% Number of intervals
lenwin = size(win,2);

% Interval start and end points
allc = [ones(1,lenwin); -1*ones(1,lenwin)];  % +1 for start and -1 for end points

% Convert to disjunct intervals
alliml = win(:);
allcl = allc(:);
ssegs = sortrows([alliml allcl],1);   % all start and end points
msegs = [ssegs(:,1) cumsum(ssegs(:,2))];   % non-overlapping intervals start at +1 (that come after a 0) and end at 0
ndinx = find(msegs(:,2)==0);   % new end points
stinx = [1; ndinx(1:end-1)+1];   % new start points (first point is start; after every end point, a start point must follow (except the last one)

% Cumulative length of all intervals
len = sum(ssegs(ndinx,1)-ssegs(stinx,1));