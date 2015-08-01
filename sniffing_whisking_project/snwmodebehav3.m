function snwmodebehav3
%SNWMODEBEHAV3   Coupling mode statistics.
%   SNWMODEBEHAV3 quantifies the clustering of sniffing-whisking coupling
%   modes around behavioral events (cue tone onset). It evaluates
%   statistical significance from the point of view of the tones, i.e. is a
%   certain mode more prevalent in a window after tone onset then before.
%   The overlap between segments of a given mode (see SNWSNIFFPLOT2) with
%   windows relative to behavioral time stamps is assessed. The number of
%   behav. windows that overlap with segments of a given mode mode is
%   computed. Chi-square test is used for statistical comparison.
%
%   See also SNWMODEBEHAV, SNWSNIFFPLOT2, SNWMODESWITCH and SNWMODESWITCH2.

% Sampling rate
sr = 1000;

% Load all segments
global DATADIR
fn = [DATADIR 'HSW\allsegments_behav.mat'];
load(fn)
fn = [DATADIR 'HSW\tones.mat'];
load(fn)
fn = [DATADIR 'HSW\switches.mat'];
load(fn)

% Mode indices
gr = [allsegments.group];   % group codes corresponding to all segments
gr05 = find(gr==0.5);   % indices for 1:2 mode
gr1 = find(gr==1);   % indices for 1:1 mode
gr2 = find(gr==2);   % indices for 2:1 mode

% Segment information
rrat = {allsegments.rat};   % rat ID
rdate = {allsegments.session};   % session date

% Trial information
isvalid1 = [tones.lopresniff];
isvalid2 = [tones.hipresniff];

% Quantify grouping of modes after tone onset
NumTones = length(tones);   % aggregate number of tone presentations
[aspc05 aspc1 aspc2] = deal(nan(1,NumTones));   % preallocate
for iS = 1:NumTones   % loop through all cue tone presentations
    sessinx = find(strcmp(tones(iS).session,rdate));   % find segments from the same session as the cue tone
    seg05 = allsegments(intersect(sessinx,gr05));   % 1:2 segments from the session
    seg1 = allsegments(intersect(sessinx,gr1));   % 1:1 segments from the session
    seg2 = allsegments(intersect(sessinx,gr2));   % 2:1 segments from the session
    segwin05 = [seg05.start; seg05.end];   % first row: start points; second row: end points of the 1:2 segments
    segwin1 = [seg1.start; seg1.end];   % first row: start points; second row: end points of the 1:1 segments
    segwin2 = [seg2.start; seg2.end];   % first row: start points; second row: end points of the 2:1 segments
    tonewin = [tones(iS).timestamp+2; tones(iS).timestamp+8];   % window around tone onset
    if ~isempty(segwin05)
        aspc05(iS) = winoverlap(tonewin,segwin05);   % percent of the 1:2 segment in the tone window
    end
    if ~isempty(segwin1)
        aspc1(iS) = winoverlap(tonewin,segwin1);   % percent of the 1:1 segment in the tone window
    end
    if ~isempty(segwin2)
        aspc2(iS) = winoverlap(tonewin,segwin2);   % percent of the 2:1 segment in the tone window
    end
end

% Quantify grouping of modes before tone onset
[bspc05 bspc1 bspc2] = deal(nan(1,NumTones));   % preallocate
for iS = 1:NumTones   % loop through all cue tone presentations
    sessinx = find(strcmp(tones(iS).session,rdate));   % find segments from the same session as the cue tone
    seg05 = allsegments(intersect(sessinx,gr05));   % 1:2 segments from the session
    seg1 = allsegments(intersect(sessinx,gr1));   % 1:1 segments from the session
    seg2 = allsegments(intersect(sessinx,gr2));   % 2:1 segments from the session
    segwin05 = [seg05.start; seg05.end];   % first row: start points; second row: end points of the 1:2 segments
    segwin1 = [seg1.start; seg1.end];   % first row: start points; second row: end points of the 1:1 segments
    segwin2 = [seg2.start; seg2.end];   % first row: start points; second row: end points of the 2:1 segments
    tonewin = [tones(iS).timestamp-8; tones(iS).timestamp-2];   % window around tone onset
    if ~isempty(segwin05)
        bspc05(iS) = winoverlap(tonewin,segwin05);   % percent of the 1:2 segment in the tone window
    end
    if ~isempty(segwin1)
        bspc1(iS) = winoverlap(tonewin,segwin1);   % percent of the 1:1 segment in the tone window
    end
    if ~isempty(segwin2)
        bspc2(iS) = winoverlap(tonewin,segwin2);   % percent of the 2:1 segment in the tone window
    end
end

% Statistics
signrank(bspc05,aspc05)   % Wilcoxon sign rank tests
signrank(bspc1,aspc1)
signrank(bspc2,aspc2)

[h p] = b_chi2test2([sum(bspc05>0) sum(bspc05<=0)],[sum(aspc05>0) sum(aspc05<=0)],0.05)   % chi square test
[h p] = b_chi2test2([sum(bspc1>0) sum(bspc1<=0)],[sum(aspc1>0) sum(aspc1<=0)],0.05)
[h p] = b_chi2test2([sum(bspc2>0) sum(bspc2<=0)],[sum(aspc2>0) sum(aspc2<=0)],0.05)

% Statistics for low pre-pellet sniff frequency
signrank(bspc05(isvalid1),aspc05(isvalid1))   % restricted to valid trials
signrank(bspc1(isvalid1),aspc1(isvalid1))
signrank(bspc2(isvalid1),aspc2(isvalid1))

[h p] = b_chi2test2([sum(bspc05(isvalid1)>0) sum(bspc05(isvalid1)<=0)],[sum(aspc05(isvalid1)>0) sum(aspc05(isvalid1)<=0)],0.05)   % restricted to valid trials
[h p] = b_chi2test2([sum(bspc1(isvalid1)>0) sum(bspc1(isvalid1)<=0)],[sum(aspc1(isvalid1)>0) sum(aspc1(isvalid1)<=0)],0.05)
[h p] = b_chi2test2([sum(bspc2(isvalid1)>0) sum(bspc2(isvalid1)<=0)],[sum(aspc2(isvalid1)>0) sum(aspc2(isvalid1)<=0)],0.05)

% Statistics for high pre-pellet sniff frequency
signrank(bspc05(isvalid2),aspc05(isvalid2))   % restricted to valid trials
signrank(bspc1(isvalid2),aspc1(isvalid2))
signrank(bspc2(isvalid2),aspc2(isvalid2))

[h p] = b_chi2test2([sum(bspc05(isvalid2)>0) sum(bspc05(isvalid2)<=0)],[sum(aspc05(isvalid2)>0) sum(aspc05(isvalid2)<=0)],0.05)   % restricted to valid trials
[h p] = b_chi2test2([sum(bspc1(isvalid2)>0) sum(bspc1(isvalid2)<=0)],[sum(aspc1(isvalid2)>0) sum(aspc1(isvalid2)<=0)],0.05)
[h p] = b_chi2test2([sum(bspc2(isvalid2)>0) sum(bspc2(isvalid2)<=0)],[sum(aspc2(isvalid2)>0) sum(aspc2(isvalid2)<=0)],0.05)

keyboard

% Statistics comparing low and high pre-pellet frequencies
[h p] = b_chi2test2([sum(bspc05(isvalid1)>0) sum(bspc05(isvalid1)<=0)],[sum(bspc05(isvalid2)>0) sum(bspc05(isvalid2)<=0)],0.05)
[h p] = b_chi2test2([sum(bspc1(isvalid1)>0) sum(bspc1(isvalid1)<=0)],[sum(bspc1(isvalid2)>0) sum(bspc1(isvalid2)<=0)],0.05)
[h p] = b_chi2test2([sum(bspc2(isvalid1)>0) sum(bspc2(isvalid1)<=0)],[sum(bspc2(isvalid2)>0) sum(bspc2(isvalid2)<=0)],0.05)

% -------------------------------------------------------------------------
function [spc ppc] = winoverlap(segwin,tonewin)

% Calculate length covered by the segments
allwin = [tonewin segwin];   % "union" of segments
lpw = lebesgue(tonewin);  % tone windows
lsw = lebesgue(segwin);  % current segment
law = lebesgue(allwin);   % all segments ("union")

% Calculate overlap
olp = lpw + lsw - law;   % calculate intersect from union

%  Relative overlap
spc = olp / lsw;   % overlap relative to the current segment
ppc = olp / lpw;   % overlap relative to the tone windows

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