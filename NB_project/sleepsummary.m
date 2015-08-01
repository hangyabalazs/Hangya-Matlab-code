function sleepsummary
%SLEEPSUMMARY   Firing rate during sleep.
%   SLEEPSUMMARY performs summary analysis for freely moving/sleep
%   recordings. Firing rate during sleep, quiet wakefulness and freely
%   moving epoch is plotted. Firing rates are compared with Mann-Whitney
%   U-test (firing rates binned in 5s segments for the test only).
%
%   See also FIRINGRATE_ANALYSIS_SLEEP.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   11-Dec-2012

%   Edit log: BH, 12/11/13

% All ChAT and pChAT cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
allChAT = [ChAT pChAT];

% n029_120203a_3.1
[FiringRate(1,1:3) FiringRateBinned(1,1:3) Speed(1,1:3)] = firingrate_analysis_sleep(allChAT{3});
FiringRate(1,4) = getvalue('baseline','n029_120203a_3.1');

% n046_130101a_6.1
[FiringRate(2,1:3) FiringRateBinned(2,1:3) Speed(2,1:3)] = firingrate_analysis_sleep(allChAT{9});
FiringRate(2,4) = getvalue('baseline','n046_130101a_6.1');

% n046_130104a_6.2
[FiringRate(3,1:3) FiringRateBinned(3,1:3) Speed(3,1:3)] = firingrate_analysis_sleep(allChAT{12});
FiringRate(3,4) = getvalue('baseline','n046_130104a_6.2');

% n046_130108a_4.1
[FiringRate(4,1:3) FiringRateBinned(4,1:3) Speed(4,1:3)] = firingrate_analysis_sleep(allChAT{13});
FiringRate(4,4) = getvalue('baseline','n046_130108a_4.1');

% n046_130108a_4.2
[FiringRate(5,1:3) FiringRateBinned(5,1:3) Speed(5,1:3)] = firingrate_analysis_sleep(allChAT{14});
FiringRate(5,4) = getvalue('baseline','n046_130108a_4.2');

% Plot
figure
hold on
for iC = 1:5
    for x1 = 1:2
        x2 = x1 + 1;
        p = ranksum(FiringRateBinned{iC,x1},FiringRateBinned{iC,x2});
        if p < 0.01
            line([x1 x2],[FiringRate(iC,x1) FiringRate(iC,x2)],'Color',[0.7 0.7 0.7])
        else
            line([x1 x2],[FiringRate(iC,x1) FiringRate(iC,x2)],...
                'Color',[0.7 0.7 0.7],'LineStyle',':')
        end
    end
end
bar(1:3,median(FiringRate(:,1:3)),'BarWidth',0.5,'EdgeColor','k',...
    'FaceColor','none','LineWidth',2)
keyboard

% -------------------------------------------------------------------------
function [FiringRate FiringRateBinned Speed] = firingrate_analysis_sleep(cellid,varargin)
%FIRINGRATE_ANALYSIS_SLEEP   Firing rate.
%   FR = FIRINGRATE_ANALYSIS_SLEEP(CELLID) calculates firing rate in sleep 
%   sessions. LFP spectrum, delta to all power ratio, temporal evolution of
%   firing rate and the animals speed, speed - firing rate association are
%   plotted.
%
%   Note that the sleep session MAY INCLUDE LIGHT_STIMULATION EPISODES!
%
%   See also FIRINGRATE_ANALYSIS and NBFIRINGRATE.

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
g = parse(prs,cellid,varargin{:});

% Load sleep segments
global DATAPATH
inpdir = fullfile(DATAPATH,'NB','speed',filesep);
fnm = fullfile(inpdir,'sleep_segments2.xlsx');
[xnum xstr xcell] = xlsread(fnm);     % read Excel file
[rowinx, ~] = find(strcmp(cellid,xcell));   % find cell in table
SegNum = length(rowinx);   % number of segments

% Load session data
sleepsessinx = strcmp('sleep cellid',xcell(1,:));
for iS = 1:SegNum   % loop through segments
    sleepcellid = xcell{rowinx(iS),sleepsessinx};   % sleep session
    sleepcellidt = regexprep(sleepcellid,'\.','_');
    if ~exist('SP','var') || ~isfield(SP,sleepcellidt)
        [SP.(sleepcellidt),spd.(sleepcellidt),time_speed.(sleepcellidt),stim_segs.(sleepcellidt)] = ...
            get_session_data(sleepcellid);  % load spike and position data, return speed
    end
end

% Freely moving
[FiringRate(1) FiringRateBinned{1} Speed(1)] = segment_speed(cellid,SP,spd,time_speed,stim_segs,xcell,'freely moving');

% Quiet wakefulness
[FiringRate(2) FiringRateBinned{2} Speed(2)] = segment_speed(cellid,SP,spd,time_speed,stim_segs,xcell,'quiet wakefulness');

% Sleep
[FiringRate(3) FiringRateBinned{3} Speed(3)] = segment_speed(cellid,SP,spd,time_speed,stim_segs,xcell,'sleep');
 
% -------------------------------------------------------------------------
function [FiringRate FiringRateBinned Speed] = segment_speed(cellid,SP,spd,time_speed,stim_segs,xcell,str)

% Firing rate
colinx = find(strcmp(str,xcell(1,:)));   % find segment boundaries in table
[rowinx, ~] = find(strcmp(cellid,xcell));   % find cell in table
SegNum = length(rowinx);   % number of segments[pfr psp] = deal(nan(1,SegNum));
[pfr psp] = deal(nan(1,SegNum));
FiringRateBinned = [];
sleepsessinx = strcmp('sleep cellid',xcell(1,:));
for iS = 1:SegNum   % loop through segments
    sleepcellid = xcell{rowinx(iS),sleepsessinx};   % cell ID for the sleep session
    sleepcellidt = regexprep(sleepcellid,'\.','_');
    lSP = SP.(sleepcellidt);   % spike data of the sleep session
    lspd = spd.(sleepcellidt);   % speed data of the sleep session
    ltime_speed = time_speed.(sleepcellidt);   % time vector for speed data
    lstim_segs = stim_segs.(sleepcellidt);   % stimulation segments
    start_time = xcell{rowinx(iS),colinx};   % segments start time(s)
    end_time = xcell{rowinx(iS),colinx+1};   % segment end time(s)
    tseg = [start_time lstim_segs(2,:); lstim_segs(1,:) end_time];   % gaps between stim. segments
    if any(lstim_segs(:)>start_time&lstim_segs(:)<end_time)
        error('Stimulation has to be excluded.')
    end
    
    pfr(iS) = sum(lSP>start_time&lSP<end_time) / ...
        (end_time - start_time);   % firing rate for the segment
    psp(iS) = mean(lspd(find(ltime_speed>start_time,1,'first'):...
        find(ltime_speed<end_time,1,'last')));   % speed for the segment
    
    % Firing rate binned to 5s segments
    if ~isnan(start_time)
        seglen = 5;        % 5 sec. long segments
        len = end_time - start_time;   % length of full epoch
        lim1 = start_time:seglen:end_time-seglen;   % segment boundaries
        lim2 = lim1 + seglen;
        SegNum2 = floor(len/seglen);   % number of segments
        pfr2 = nan(1,SegNum2);
        for iS2 = 1:SegNum2
            segspikes = lSP(lSP>lim1(iS2)&lSP<lim2(iS2));  % spikes in the segment
            pfr2(iS2) = length(segspikes) / seglen;  % FR in the segment
        end
        FiringRateBinned = [FiringRateBinned pfr2];  % binned FR
    end
end
FiringRate = nanmean(pfr);   % average FR over segments
Speed = nanmean(psp);   % average speed over segments

% -------------------------------------------------------------------------
function [SP,spd,time_speed,stim_segs] = get_session_data(cellid)

% Load absoulte spike times
SP = loadcb(cellid);   % absoulte spike times

% Position data
[animalID sessionID] = cellid2tags(cellid);
fullpth = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
fn = fullfile(fullpth,'VT1.nvt');
[TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points, NlxHeader] =...
    Nlx2MatVT(fn,[1 1 1 1 1 1],1,1);   % load position data
TimeStamps = TimeStamps / 1e6;   % convert time stamps to seconds

% Position correction
[xpos ypos time_pos] = poscorrection(ExtractedX,ExtractedY,TimeStamps);
if ~isequal(time_pos,TimeStamps)  % internal check
    error('Bug.')
end

% Resample position vector
dt_pos = mean(diff(time_pos));
sr_pos = 1 / dt_pos;    % original sampling rate (changed later)

% Speed
rx = diff(prctile(xpos,[0.25 99.75]));   % x range, allowing some outliers
ry = diff(prctile(ypos,[0.25 99.75]));   % y range, allowing some outliers
xpos2 = xpos / rx * 27.4;  % x position in cm
ypos2 = ypos / ry * 16.5;  % y position in cm
spd = mouse_speed(xpos2,ypos2,sr_pos);   % speed vector
time_speed = time_pos(1:end-1);   % corresponding time vector

% Load Events
% fn = [fullpth '\EVENTS.mat'];
% load(fn);         % load converted Neuralynx events
% figure;
% plot(Events_TimeStamps,Events_Nttls)

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');
    
    % Druing stimulation protocols: from pulse onset
    stim_segs = [SE.PulseOn - 1; SE.PulseOn + 5];
catch
    disp([cellid ': There was no stim protocol for this session.'])
    stim_segs = [];
end

% -------------------------------------------------------------------------
function [X Y TP] = poscorrection(x_pos,y_pos,x_pos_ts)

% Position vector correction
bns = 8;   % spatial bin size
TP2 = x_pos_ts;   % for later interpolation
dx = diff(x_pos);       % leave out the outliers
fdx = find(abs(dx)>bns&abs([dx(2:end) 0])>bns&abs(dx-[dx(2:end) 0])>2*bns);
x_pos(fdx+1) = [];
y_pos(fdx+1) = [];
x_pos_ts(fdx+1) = [];
dy = diff(y_pos);
fdy = find(abs(dy)>bns&abs([dy(2:end) 0])>bns&abs(dy-[dy(2:end) 0])>2*bns);
x_pos(fdy+1) = [];
y_pos(fdy+1) = [];
x_pos_ts(fdy+1) = [];

inxs = x_pos > 0 & y_pos > 0;   % leave out (0;0) points
X = x_pos(inxs);
Y = y_pos(inxs);
TP = x_pos_ts(inxs);

% Interpolate the deleted outlier timepoints.
X = interp1(TP,X,TP2);
Y = interp1(TP,Y,TP2);
TP = TP2;

% -------------------------------------------------------------------------
function spd = mouse_speed(xpos,ypos,sr)

% Speed
scxp = smooth(xpos,'linear',51);   % smoothed x positions
scyp = smooth(ypos,'linear',51);   % smoothed y positions
dx = diff(scxp);   % derivative of x positions
dy = diff(scyp);   % derivative of y positions
spd = sqrt(dx.^2+dy.^2) * sr;   % speed