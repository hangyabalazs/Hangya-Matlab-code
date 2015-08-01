function FAlatjit_HDB
%FALATJIT   Latency and jitter.
%   FALATJIT calculates spike latency and jitter distribution of
%   cholinergic neurons with respect to punishment. Latency is calculated
%   as PSTH peak and jitter is computed based on first spike times within a
%   detected response window (see NBACTIVATION for details). The cumulative
%   density functions are plotted and saved.
%
%   See also NBACTIVATION and HIT_VS_FARESPONSE.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   2-Oct-2013

%   Edit log: BH 8/2/13

% Directories
global DATAPATH
resdir = [DATAPATH 'HDB\FAlatjit\ChAT_newdata\'];
issave = true;

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified (n = 12)
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
allChAT = ChAT;   % identified and putative

% Response properties
FA_latency = getvalue('FA_Latency',allChAT);   % latency for false alarms
Hit_latency = getvalue('Hit_Latency',allChAT);   % latency for hits
FA_jitter = getvalue('FA_Jitter',allChAT);   % jitter for hits
Hit_jitter = getvalue('Hit_Jitter',allChAT);   % jitter for false alarms
FA_reliability = getvalue('FA_Reliability',allChAT);   % reliability for false alarms
Hit_reliability = getvalue('Hit_Reliability',allChAT);   % reliability for hits
FA_spno = getvalue('FA_SpikeNumberDistribution',allChAT);   % intra-burst spike number for false alarms
Hit_spno = getvalue('Hit_SpikeNumberDistribution',allChAT);   % intra-burst spike number for hits
FA_spno = cellfun(@(s)mean(s(s>0)),FA_spno);   % average spike number
Hit_spno = cellfun(@(s)mean(s(s>0)),Hit_spno);

% PSTH statistics
FA_psth_stats = getvalue('FA_psth_stats',allChAT);
Hit_psth_stats = getvalue('Hit_psth_stats',allChAT);
FA_psth_stats = nancell2struct(FA_psth_stats);
Hit_psth_stats = nancell2struct(Hit_psth_stats);
FA_p = round([FA_psth_stats.Wpa].*1000) / 1000;   % round p value to 1/1000
Hit_p = round([Hit_psth_stats.Wpa].*1000) / 1000;   % round p value to 1/1000
FA_act = FA_p <= 0.025;  % cells activated after false alarms
Hit_act = Hit_p <= 0.025;   % cells activated after hits

% FA latency CDF
inx = FA_act;  % cells activated in false alarm trials
red = [216 41 0] / 255;
figure
edges = 0:0.0001:max(FA_latency(inx))+0.01;   % bin edges
dist_latency = histc(FA_latency(inx),edges);   % latency histogram
dist_latency = [0; dist_latency(1:end-1)];   % values corresponding to the edges
dist_latency = dist_latency / sum(dist_latency);   % normalize
stairs(edges,cumsum(dist_latency),'Color','red','LineWidth',2)
axis tight
set(gca,'box','off','FontSize',12,'TickDir','out')
xlabel('Latency')
setmyplot_balazs

fnm = [resdir 'FA_latency_CDF.fig'];   % save
if issave
    saveas(gcf,fnm)
end

% FA jitter CDF
inx = FA_act;  % cells activated in false alarm trials
figure
edges = 0:0.0001:max(FA_jitter(inx)+0.001);   % bin edges
dist_jitter = histc(FA_jitter(inx),edges);   % jitter histogram
dist_jitter = [0; dist_jitter(1:end-1)];   % values corresponding to the edges
dist_jitter = dist_jitter / sum(dist_jitter);   % normalize
stairs(edges,cumsum(dist_jitter),'Color','red','LineWidth',2)
axis tight
set(gca,'box','off','FontSize',12,'TickDir','out')
xlabel('Jitter')
setmyplot_balazs

fnm = [resdir 'FA_jitter_CDF.fig'];   % save
if issave
    saveas(gcf,fnm)
end

% FA reliability CDF
inx = FA_act;  % cells activated in false alarm trials
figure
edges = 0:0.0001:max(FA_reliability(inx)+0.001);   % bin edges
dist_reliability = histc(FA_reliability(inx),edges);   % reliability histogram
dist_reliability = [0; dist_reliability(1:end-1)];   % values corresponding to the edges
dist_reliability = dist_reliability / sum(dist_reliability);   % normalize
stairs(edges,cumsum(dist_reliability),'Color','red','LineWidth',2)
axis tight
set(gca,'box','off','FontSize',12,'TickDir','out')
xlabel('Reliability')
setmyplot_balazs

fnm = [resdir 'FA_reliability_CDF.fig'];   % save
if issave
    saveas(gcf,fnm)
end

% FA intra-burst spike number CDF
inx = FA_act;  % cells activated in false alarm trials
figure
edges = 0:0.0001:max(FA_spno(inx)+0.001);   % bin edges
dist_spno = histc(FA_spno(inx),edges);   % intra-burst spike number histogram
dist_spno = [0; dist_spno(1:end-1)];   % values corresponding to the edges
dist_spno = dist_spno / sum(dist_spno);   % normalize
stairs(edges,cumsum(dist_spno),'Color','red','LineWidth',2)
axis tight
set(gca,'box','off','FontSize',12,'TickDir','out')
xlabel('Number of spikes in burst')
setmyplot_balazs

fnm = [resdir 'FA_spno_CDF.fig'];   % save
if issave
    saveas(gcf,fnm)
end