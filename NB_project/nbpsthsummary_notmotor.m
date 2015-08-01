function nbpsthsummary_notmotor
%NBPSTHSUMMARY_NOTMOTOR   PSTH summary figures.
%   NBPSTHSUMMARY_NOTMOTOR plots population PSTH for cholinergic neurons
%   with feedback-delay. Z-scored PSTHs are used throughout the analysis.
%   Average PSTHs aligned to reponse and feedback are plotted after
%   baseline-subtraction.
%
%   See also NBPSTHSUMMARY.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   8-March-2014

%   Edit log: BH 3/8/14

% Directories
global DATAPATH
resdir = ([DATAPATH 'NB\psthsummary_notmotor_newdata\']);   % result directory

% Cells
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % ChAT cells with FeedbackDelay
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
NumCells = length(ChAT);

% Load data
allpsth = getvalue('FA_psth',ChAT);
allpsth = nancell2mat(allpsth);
allpsth = zscore(allpsth,0,2);
allstats = getvalue('FA_psth_stats',ChAT);
allstats = nancell2struct2(allstats);

% LeftPortIn PSTH
allpsth_lpi = [];
wn = [-0.6 0.6];
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector
for k = 1:NumCells
    cellid = ChAT{k};   % cell ID
    disp(cellid)
    
    % Align event
    alignfilter = 'FalseAlarm==1';
    alignevent = 'LeftPortIn';
    
    % Calcualte PSTH
    [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
        alignevent,wn,...
        'dt',dt,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
        'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
    
    % Concatenate data from different cells
    allpsth_lpi = [allpsth_lpi; psth]; %#ok<AGROW>
end

% Normalize
allpsth = zscore(allpsth,0,2);
allpsth_lpi = zscore(allpsth_lpi,0,2);

% Population PSTHs for groups
figure   % images proportional to number of cells in the groups
imagesc(allpsth);
colormap hot
set(gca,'CLim',[-2 20])
saveas(gcf,fullfile(resdir,'poppsth_ChAT_fb.fig'))

figure
imagesc(allpsth_lpi);
colormap hot
set(gca,'CLim',[-2 20])
saveas(gcf,fullfile(resdir,'poppsth_ChAT_lpi.fig'))

% Average PSTH
figure
hold on;
baseline_ChAT_fb = mean(allpsth(:,1:200));
mn_fb = mean(baseline_ChAT_fb);   % baseline for feedback alignment
errorshade(time,mean(allpsth)-mn_fb,std(allpsth)/sqrt(size(allpsth,1)),...
    'LineColor','k','ShadeColor','k')

ymx = 12;
set(gca,'YLim',[-1 ymx],'YTick',[0 ymx/2 ymx],'YTickLabel',{'0' '' num2str(ymx)},...
    'XLim',[-0.6 0.6],'XTick',[-0.6 -0.3 0 0.3 0.6],'XTickLabel',{'-600' '' '0' '' '600'});
line([0 0],ylim,'Color','k','LineStyle',':')
line(xlim,[0 0],'Color','k')
xlabel('Time from punishment (ms)')
ylabel({'Normalized';'firing rate'})
setmyplot_balazs
saveas(gcf,fullfile(resdir,'psth_average_fb.fig'))
set(gcf,'Renderer','painters')
saveas(gcf,fullfile(resdir,'psth_average_fb.eps'))

figure
hold on;
baseline_ChAT_lpi = mean(allpsth_lpi(:,1:200));
mn_lpi = mean(baseline_ChAT_lpi);   % baseline for response alignment
errorshade(time,mean(allpsth_lpi)-mn_lpi,std(allpsth_lpi)/sqrt(size(allpsth_lpi,1)),...
    'LineColor','k','ShadeColor','k')

set(gca,'YLim',[-1 ymx],'YTick',[0 ymx/2 ymx],'YTickLabel',{'0' '' num2str(ymx)},...
    'XLim',[-0.6 0.6],'XTick',[-0.6 -0.3 0 0.3 0.6],'XTickLabel',{'-600' '' '0' '' '600'});
line([0 0],ylim,'Color','k','LineStyle',':')
line(xlim,[0 0],'Color','k')
xlabel('Time from response (ms)')
ylabel({'Normalized';'firing rate'})
setmyplot_balazs
saveas(gcf,fullfile(resdir,'psth_average_lpi.fig'))
set(gcf,'Renderer','painters')
saveas(gcf,fullfile(resdir,'psth_average_lpi.eps'))

keyboard