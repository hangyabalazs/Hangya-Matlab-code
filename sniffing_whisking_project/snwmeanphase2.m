function snwmeanphase2
%SNWMENPHASE2   Whisking phase relative to respiration.
%   SNWMENPHASE2 calculates summary data from the phase statistics generated
%   by SNWSNIFFPLOT (see SNWSNIFFPLOT for details). Mean whisking
%   conditioned on sniffing phase is plotted for different coulping types.
%   Mean is calculated for individual animals first.
%
%   See also SNWMEANPHASE and SNWSNIFFPLOT.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Phase
rats = unique(tbl(:,1));
numRats = length(rats);
allang1to1 = [];
allang1to2 = [];
allang1to2b = [];
allang2to1 = [];
allang1to1split = [];
allang1to1splitb = [];
for iR = 1:numRats
    crat = rats(iR);
    sessinx = find(strcmp(tbl(:,1),crat));
    
    numSessions = length(sessinx);
    callang1to1 = [];
    callang1to2 = [];
    callang1to2b = [];
    callang2to1 = [];
    callang1to1split = [];
    callang1to1splitb = [];
    for iS = 1:numSessions
        rrat = tbl{sessinx(iS),1};
        rdate = tbl{sessinx(iS),2};
        [cnts aang1to1 aang1to2 aang1to2b aang2to1 aang1to1split aang1to1splitb]...
            = main(rrat,rdate);
        close all
        callang1to1 = [callang1to1; aang1to1/sum(aang1to1)];
        callang1to2 = [callang1to2; aang1to2/sum(aang1to2)];
        callang1to2b = [callang1to2b; aang1to2b/sum(aang1to2)];
        callang2to1 = [callang2to1; aang2to1/sum(aang2to1)];
        callang1to1split = [callang1to1split; aang1to1split/sum(aang1to1split)];
        callang1to1splitb = [callang1to1splitb; aang1to1splitb/sum(aang1to1split)];
    end
    allang1to1 = [allang1to1; nanmean(callang1to1,1)];
    allang1to2 = [allang1to2; nanmean(callang1to2,1)];
    allang1to2b = [allang1to2b; nanmean(callang1to2b,1)];
    allang2to1 = [allang2to1; nanmean(callang2to1,1)];
    allang1to1split = [allang1to1split; nanmean(callang1to1split,1)];
    allang1to1splitb = [allang1to1splitb; nanmean(callang1to1splitb,1)];
end
keyboard

figure
plot([cnts cnts+2*pi],[nanmean(allang1to1) nanmean(allang1to1)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang1to1) / sqrt(size(allang1to1,1));
errorshade([cnts cnts+2*pi],[nanmean(allang1to1) nanmean(allang1to1)],[tse tse],'r-')

figure
plot([cnts cnts+2*pi],[nanmean(allang1to1split) nanmean(allang1to1split)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang1to1split) / sqrt(size(allang1to1split,1));
errorshade([cnts cnts+2*pi],[nanmean(allang1to1split) nanmean(allang1to1split)],[tse tse],...
    'LineColor','r','ShadeColor','r')
tse = nanstd(allang1to1splitb) / sqrt(size(allang1to1splitb,1));
errorshade([cnts cnts+2*pi],[nanmean(allang1to1splitb) nanmean(allang1to1splitb)],[tse tse],...
    'LineColor','b','ShadeColor','b')

figure
plot([cnts cnts+2*pi],[nanmean(allang1to2) nanmean(allang1to2)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang1to2) / sqrt(size(allang1to2,1));
errorshade([cnts cnts+2*pi],[nanmean(allang1to2) nanmean(allang1to2)],[tse tse],...
    'LineColor','r','ShadeColor','r')
tse = nanstd(allang1to2b) / sqrt(size(allang1to2b,1));
errorshade([cnts cnts+2*pi],[nanmean(allang1to2b) nanmean(allang1to2b)],[tse tse],...
    'LineColor','b','ShadeColor','b')

figure
plot([cnts cnts+2*pi],[nanmean(allang2to1) nanmean(allang2to1)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang2to1) / sqrt(size(allang2to1,1));
errorshade([cnts cnts+2*pi],[nanmean(allang2to1) nanmean(allang2to1)],[tse tse],'r-')

figure
plot([cnts cnts+2*pi],[smooth(nanmean(allang2to1),'linear',3) smooth(nanmean(allang2to1),'linear',3)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang2to1) / sqrt(size(allang2to1,1));
errorshade([cnts cnts+2*pi],[smooth(nanmean(allang2to1),'linear',3) smooth(nanmean(allang2to1),'linear',3)],...
    [smooth(tse,'linear',3) smooth(tse,'linear',3)],'r-')

% -------------------------------------------------------------------------
function [cnts aang1to1 aang1to2 aang1to2b aang2to1 aang1to1split aang1to1splitb]...
    = main(rrat,rdate)

% Define results directory
global DATAPATH
global DATADIR
resdir = [DATAPATH 'SniffWhisk\phasehist_summary\'];

% Load data for phase calculation
fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE_' rrat '_' rdate '.mat'];
load(fnme)
fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE1TO1_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
[R p cnts aang1to1] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase

% Load data for phase calculation
fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE1TO2_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
[R p cnts aang1to2] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
[R p cnts aang1to2b] = dphist(time_grp2,rms_whisk,resp_phase2);

% Load data for phase calculation
fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE2TO1_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
[R p cnts aang2to1] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase

% Load data for phase calculation
fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE1TO1SPLIT_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
[R p cnts aang1to1split] = dphist(time_grp_split1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
[R p cnts aang1to1splitb] = dphist(time_grp_split2,rms_whisk,resp_phase2);

% keyboard

% -------------------------------------------------------------------------
function [R p cnts mn_wh] = dphist(time_grp1,rms_whisk,resp_phase2) 
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = rms_whisk' .* time_grp1;
zeroinx = time_grp1==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
if ~isempty(find(time_grp1,1))
    [R p] = lincirc_corr2(mn_wh',cnts');
else
    R = NaN;
    p = NaN;
end