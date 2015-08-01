function nbpsthsummary_HDB
%NBPSTHSUMMARY   PSTH summary figures.
%   NBPSTHSUMMARY plots population PSTH for activated cholinergic
%   (including putative), activated non-tagged and inhibited non-tagged
%   cells. Z-scored PSTHs are used throughout the analysis. Avergae PSTHs
%   are also plotted after baseline-subtraction.
%
%   Next, average PSTH aligned to reward or punishment is plotted for all
%   (not only significantly modulated) ChAT+ cells. A heat map of the
%   population PSTH for the significantly activated (p<0.05, two-sided
%   Mann-Whitney test) ChAT+ cells is also plotted separately for hit and
%   false alarm trials.
%
%   See also PLOT_POPSTH_CELL_NEW.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   29-Sept-2013

%   Edit log: BH 9/29/13

% Directories
global DATAPATH
resdir = ([DATAPATH 'HDB\psthsummary_newdata\']);   % result directory

% Cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
NumChAT = length(ChAT);   % number of cholinergic cells

% Time
wn = [-600 600];
dt = 1;
time = wn(1):dt:wn(2);   % time vector

% Load data
allpsth = getvalue('FA_psth',ChAT);
allpsth = nancell2mat(allpsth);
FA_allpsth = zscore(allpsth,0,2);
allstats = getvalue('FA_psth_stats',ChAT);
allstats = nancell2struct2(allstats);
psthplot(resdir,time,FA_allpsth,'hit')

% Load data
allpsth = getvalue('Hit_psth',ChAT);
allpsth = nancell2mat(allpsth);
Hit_allpsth = zscore(allpsth,0,2);
allstats = getvalue('Hit_psth_stats',ChAT);
allstats = nancell2struct2(allstats);
psthplot(resdir,time,Hit_allpsth,'hit')

% Sorted popPSTH
[m1 m2] = max(FA_allpsth,[],2);
[srt Ia] = sort(m2,'ascend');   % sort based on FA response
figure   % plot all FA PSTHs, sorted
imagesc(time,1:size(FA_allpsth,1),FA_allpsth(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'ChAT_poppsth_FA_all_sorted_new.fig'))

figure   % plot all Hit PSTHs; sort based on FA response
imagesc(time,1:size(Hit_allpsth,1),Hit_allpsth(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'ChAT_poppsth_Hit_all_sorted_new.fig'))

% Average PSTH
figure   % plot
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
baseline_FA = mean(FA_allpsth(:,1:200));
mn_FA = mean(baseline_FA);   % baseline for feedback alignment
baseline_Hit = mean(Hit_allpsth(:,1:200));
mn_Hit = mean(baseline_Hit);   % baseline for feedback alignment
errorshade(time,mean(FA_allpsth)-mn_FA,std(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
hold on
errorshade(time,mean(Hit_allpsth)-mn_Hit,std(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
saveas(gcf,fullfile(resdir,['average_PSTH.fig']))

% -------------------------------------------------------------------------
function psthplot(resdir,time,allpsth,tag)

% Normalize
allpsth = zscore(allpsth,0,2);

% Population PSTHs for groups
figure   % images proportional to number of cells in the groups
imagesc(time,1:size(allpsth,1),allpsth);
colormap hot
set(gca,'CLim',[-2 20])
saveas(gcf,fullfile(resdir,['poppsth_ChAT_' tag '.fig']))

% Average PSTH
% figure
% hold on;
% baseline_ChAT = mean(allpsth(:,1:200));
% mn = mean(baseline_ChAT);   % baseline for feedback alignment
% errorshade(time,mean(allpsth)-mn,std(allpsth)/sqrt(size(allpsth,1)),...
%     'LineColor','k','ShadeColor','k')
% 
% ymx = 14;
% set(gca,'YLim',[-1 ymx],'YTick',[0 ymx/2 ymx],'YTickLabel',{'0' '' num2str(ymx)},...
%     'XLim',[-0.6 0.6],'XTick',[-0.6 -0.3 0 0.3 0.6],'XTickLabel',{'-600' '' '0' '' '600'});
% line([0 0],ylim,'Color','k','LineStyle',':')
% line(xlim,[0 0],'Color','k')
% xlabel('Time from punishment (ms)')
% ylabel({'Normalized';'firing rate'})
% setmyplot_balazs
% saveas(gcf,fullfile(resdir,['psth_average_' tag '.fig']))
% saveas(gcf,fullfile(resdir,['psth_average_' tag '.eps']))