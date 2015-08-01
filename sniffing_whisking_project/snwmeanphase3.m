function snwmeanphase3
%SNWMEANPHASE3   Whisking phase relative to respiration.
%   SNWMEANPHASE3 calculates summary data from the phase statistics generated
%   by SNWSNIFFPLOT (see SNWSNIFFPLOT for details). Mean whisking
%   conditioned on sniffing phase is plotted for different coulping types.
%
%   SNWMEANPHASE3 pools sessions from individual animals and calculates
%   mean over animals. Histograms are plotted with overlapping bins (60
%   degrees bin width, 50 degrees overlap).
%
%   See also SNWSNIFFPLOT.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Import
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
        [aang1to1 aang1to2 aang1to2b aang2to1 aang1to1split aang1to1splitb]...
            = data_import(rrat,rdate);
        close all
        callang1to1 = [callang1to1; aang1to1];
        callang1to2 = [callang1to2; aang1to2];
        callang1to2b = [callang1to2b; aang1to2b];
        callang2to1 = [callang2to1; aang2to1];
        callang1to1split = [callang1to1split; aang1to1split];
        callang1to1splitb = [callang1to1splitb; aang1to1splitb];
    end
    
    % Respiration phase - whisking relationship
    [R p cnts aang1to1] = dphist(callang1to1(:,1),callang1to1(:,2),callang1to1(:,3));    % whisking RMS conditioned on resp. phase
    
    % Respiration phase - whisking relationship
    [R p cnts aang1to2] = dphist(callang1to2(:,1),callang1to2(:,2),callang1to2(:,3));    % whisking RMS conditioned on resp. phase
    [R p cnts aang1to2b] = dphist(callang1to2b(:,1),callang1to2b(:,2),callang1to2b(:,3));
    
    % Respiration phase - whisking relationship
    [R p cnts aang2to1] = dphist(callang2to1(:,1),callang2to1(:,2),callang2to1(:,3));    % whisking RMS conditioned on resp. phase
    
    % Respiration phase - whisking relationship
    [R p cnts aang1to1split] = dphist(callang1to1split(:,1),callang1to1split(:,2),callang1to1split(:,3));    % whisking RMS conditioned on resp. phase
    [R p cnts aang1to1splitb] = dphist(callang1to1splitb(:,1),callang1to1splitb(:,2),callang1to1splitb(:,3));
    
    allang1to1 = [allang1to1; aang1to1/sum(aang1to1)];
    allang1to2 = [allang1to2; aang1to2/sum(aang1to2)];
    allang1to2b = [allang1to2b; aang1to2b/sum(aang1to2)];
    allang2to1 = [allang2to1; aang2to1/sum(aang2to1)];
    allang1to1split = [allang1to1split; aang1to1split/sum(aang1to1split)];
    allang1to1splitb = [allang1to1splitb; aang1to1splitb/sum(aang1to1split)];
    
%     keyboard

end

keyboard

% Plot
figure
plot([cnts cnts+2*pi],[nanmean(allang1to1) nanmean(allang1to1)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang1to1) / sqrt(size(allang1to1,1));
errorshade([cnts cnts+2*pi],[nanmean(allang1to1) nanmean(allang1to1)],[tse tse],...
    'LineColor',[0 153 0]/255,'ShadeColor',[0 153 0]/255)
set(gca,'FontSize',16,'TickDir','out','YLim',[0 0.06],'XLim',[cnts(1) cnts(end)+2*pi+pi/18],...
    'Box','off','XTick',[-pi 0 pi 2*pi 3*pi],'XTickLabel',[-180 0 180 360 540],...
    'YTick',[0 0.05])

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
tse = nanstd(allang1to2) / sqrt(size(allang1to2,1)-1);
errorshade([cnts cnts+2*pi],[nanmean(allang1to2) nanmean(allang1to2)],[tse tse],...
    'LineColor',[153 0 255]/255,'ShadeColor',[153 0 255]/255)
tse = nanstd(allang1to2b) / sqrt(size(allang1to2b,1)-1);
errorshade([cnts cnts+2*pi],[nanmean(allang1to2b) nanmean(allang1to2b)],[tse tse],...
    'LineColor','k','ShadeColor','k')
set(gca,'FontSize',16,'TickDir','out','YLim',[0 0.06],'XLim',[cnts(1) cnts(end)+2*pi+pi/18],...
    'Box','off','XTick',[-pi 0 pi 2*pi 3*pi],'XTickLabel',[-180 0 180 360 540],...
    'YTick',[0 0.05])

figure
plot([cnts cnts+2*pi],[nanmean(allang2to1) nanmean(allang2to1)],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
hold on
tse = nanstd(allang2to1) / sqrt(size(allang2to1,1));
errorshade([cnts cnts+2*pi],[nanmean(allang2to1) nanmean(allang2to1)],[tse tse],...
    'LineColor',[0 153 255]/255,'ShadeColor',[0 153 255]/255)
set(gca,'FontSize',16,'TickDir','out','YLim',[0 0.06],'XLim',[cnts(1) cnts(end)+2*pi+pi/18],...
    'Box','off','XTick',[-pi 0 pi 2*pi 3*pi],'XTickLabel',[-180 0 180 360 540],...
    'YTick',[0 0.05])



% -------------------------------------------------------------------------
function [aang1to1 aang1to2 aang1to2b aang2to1 aang1to1split aang1to1splitb]...
    = data_import(rrat,rdate)

% Define results directory
global DATAPATH
global DATADIR
inpdir = [DATAPATH 'SniffWhisk\phasehist_moreprecise2_midpoint\'];
resdir = [DATAPATH 'SniffWhisk\phasehist_summary\'];

% Load data for phase calculation
fnme = [inpdir 'PHASE_' rrat '_' rdate '.mat'];
load(fnme)
fnme = [inpdir 'PHASE1TO1_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
aang1to1 = [time_grp1 rms_whisk' resp_phase2];    % whisking RMS conditioned on resp. phase

% Load data for phase calculation
fnme = [inpdir 'PHASE1TO2_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
aang1to2 = [time_grp1 rms_whisk' resp_phase2];    % whisking RMS conditioned on resp. phase
aang1to2b = [time_grp2 rms_whisk' resp_phase2];

% Load data for phase calculation
fnme = [inpdir 'PHASE2TO1_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
aang2to1 = [time_grp1 rms_whisk' resp_phase2];    % whisking RMS conditioned on resp. phase

% Load data for phase calculation
fnme = [inpdir 'PHASE1TO1SPLIT_' rrat '_' rdate '.mat'];
load(fnme)

% Respiration phase - whisking relationship
aang1to1split = [time_grp_split1 rms_whisk' resp_phase2];    % whisking RMS conditioned on resp. phase
aang1to1splitb = [time_grp_split2 rms_whisk' resp_phase2];

% keyboard

% -------------------------------------------------------------------------
function [R p cnts mn_wh] = dphist(time_grp1,rms_whisk,resp_phase2) 
% edges = [-pi:pi/18:pi-pi/18; -pi+4*pi/18:pi/18:pi+3*pi/18];     % phase histogram bin limits
% cnts = (edges(1,:) + edges(2,:)) / 2;     % phase histogram bin centers
edges = [-pi-3*pi/18:pi/18:pi-4*pi/18; -pi+3*pi/18:pi/18:pi+2*pi/18];     % phase histogram bin limits
cnts = (edges(1,:) + edges(2,:)) / 2;     % phase histogram bin centers
% edges = -pi:2*pi/18:pi;     % phase histogram bin limits
% cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = rms_whisk .* time_grp1;
zeroinx = time_grp1==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep3(spvr,edges);
% [mn_phase mn_wh] = phasedep(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
if ~isempty(find(time_grp1,1))
    [R p] = lincirc_corr2(mn_wh',cnts');
else
    R = NaN;
    p = NaN;
end