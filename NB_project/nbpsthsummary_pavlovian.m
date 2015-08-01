function nbpsthsummary_pavlovian
%NBPSTHSUMMARY_PAVLOVIAN   PSTH summary figures.
%   NBPSTHSUMMARY_PAVLOVIAN plots population PSTH for cholinergic neurons
%   recorded outside behavior task. Z-scored PSTHs are used throughout the
%   analysis. Average PSTHs aligned to airpuff are plotted.
%
%   See also NBPSTHSUMMARY.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   8-March-2014

%   Edit log: BH 3/3/15

% Directories
global DATAPATH
resdir = ([DATAPATH 'NB\psthsummary_pavlovian_combined_\']);   % result directory

% Cells - NB
choosecb('NB')
selstr = ['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''pavlovian''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % ChAT cells with FeedbackDelay
% n = 1; 'n071_141210b_6.3'
NumCells = length(ChAT);

% Air-puff PSTH
allpsth = [];
allstats = [];
[Reliability Latency Jitter] = deal([]);
SpikeNumberDistribution = {};
wn = [-0.6 0.6];
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector
for k = 1:NumCells
    cellid = ChAT{k};   % cell ID
    disp(cellid)
    
    % Align event
    alignfilter = 'PTrial==1';
    alignevent = 'DeliverFeedback';
    
    % Calcualte PSTH
    [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
        alignevent,wn,...
        'dt',dt,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
        'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
        'event_filter','custom','filterinput',alignfilter,'isadaptive',2,...
        'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
        'jitterdefinition','burst','display',true);
    Reliability(end+1) = reliability;
    Latency(end+1) = latency;
    Jitter(end+1) = jitter;
    SpikeNumberDistribution{end+1} = spikenumberdistribution;
    
    % Concatenate data from different cells
    allpsth = [allpsth; psth]; %#ok<AGROW>
    allstats = [allstats; stats]; %#ok<AGROW>
end

% Cells - HDB
choosecb('HDB')
selstr = ['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''pavlovian''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})'];
ChAT = selectcell(selstr);   % ChAT cells with FeedbackDelay
ChAT = [ChAT 'h006_140409c_1.1'];
% n = 4; 'n070_150104b_5.2' 'n078_150104b_1.1' 'n078_150111b_3.1' 'h006_140409c_1.1'
NumCells = length(ChAT);
for k = 1:NumCells
    cellid = ChAT{k};   % cell ID
    disp(cellid)
    
    % Align event
    alignfilter = 'PTrial==1';
    alignevent = 'DeliverFeedback';
    
    % Calcualte PSTH
    [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
        alignevent,wn,...
        'dt',dt,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
        'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
    
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
        'event_filter','custom','filterinput',alignfilter,'isadaptive',2,...
        'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
        'jitterdefinition','burst','display',true);
    Reliability(end+1) = reliability;
    Latency(end+1) = latency;
    Jitter(end+1) = jitter;
    SpikeNumberDistribution{end+1} = spikenumberdistribution;
    
    % Concatenate data from different cells
    allpsth = [allpsth; psth]; %#ok<AGROW>
    allstats = [allstats; stats]; %#ok<AGROW>
end

% Cells - 'pavlovian' cellbase
choosecb('pavlovian')
ChAT = {'nb053_140430b_2.2'};
NumCells = length(ChAT);
for k = 1:NumCells
    cellid = ChAT{k};   % cell ID
    disp(cellid)
    
    % Align event
    alignfilter = 'PTrial==1';
    alignevent = 'DeliverFeedback';
    
    % Calcualte PSTH
    [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
        alignevent,wn,...
        'dt',dt,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
        'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
    
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
        'event_filter','custom','filterinput',alignfilter,'isadaptive',2,...
        'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
        'jitterdefinition','burst','display',true);
    Reliability(end+1) = reliability;
    Latency(end+1) = latency;
    Jitter(end+1) = jitter;
    SpikeNumberDistribution{end+1} = spikenumberdistribution;
    
    % Concatenate data from different cells
    allpsth = [allpsth; psth]; %#ok<AGROW>
    allstats = [allstats; stats]; %#ok<AGROW>
end

% Normalize
allpsth = zscore(allpsth,0,2);
allpsth = zscore(allpsth,0,2);

% Population PSTHs for groups
time = -600:600;
figure   % images proportional to number of cells in the groups
imagesc(allpsth);
colormap hot
set(gca,'CLim',[-1 20])
saveas(gcf,fullfile(resdir,'poppsth_ChAT_airpuff.fig'))

[m1 m2] = max(allpsth,[],2);
[srt Ia] = sort(m2,'ascend');   % sort based on response latency
figure   % plot all PSTHs, sorted
imagesc(time,1:size(allpsth,1),allpsth(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'poppsth_ChAT_airpuff_sorted.fig'))

% Average PSTH
figure
hold on;
baseline_ChAT_fb = mean(allpsth(:,1:200));
mn_fb = mean(baseline_ChAT_fb);   % baseline for feedback alignment
errorshade(time,mean(allpsth)-mn_fb,std(allpsth)/sqrt(size(allpsth,1)),...
    'LineColor','k','ShadeColor','k')

ymx = 10;
set(gca,'YLim',[-1 ymx],'YTick',[0 ymx/2 ymx],'YTickLabel',{'0' '' num2str(ymx)},...
    'XLim',[-600 600],'XTick',[-600 -300 0 300 600]);
line([0 0],ylim,'Color','k','LineStyle',':')
line(xlim,[0 0],'Color','k')
xlabel('Time from punishment (ms)')
ylabel({'Normalized';'firing rate'})
setmyplot_balazs
saveas(gcf,fullfile(resdir,'psth_average_airpuff.fig'))
set(gcf,'Renderer','painters')
saveas(gcf,fullfile(resdir,'psth_average_airpuff.eps'))



keyboard