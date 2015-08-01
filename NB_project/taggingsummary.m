function taggingsummary
%TAGGINGSUMMARY   Latency, jitter and reliability of tagging.
%   TAGGINGSUMMARY detects the peak activation after light stimulation and
%   determines a window 3% from baseline for activation. Peak-latency,
%   first-spike jitter and reliability (within this window) are calculated
%   and saved.
%
%   See RELIABILITY_LATENCY_JITTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   3-Oct-2013

%   Edit log: BH 10/3/13

% Directories
global DATAPATH
resdir = [DATAPATH 'NB\taggingsummary_newdata\PulseOn\'];
issave = true;

% Cholinergic neurons
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
NumChAT = length(ChAT);   % number of cholinergic cells

% H index distribution
nonChAT = selectcell(['"ChAT+"==0&"Hindex">=0.01&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})&' ...
    'ismember("RatId",[8 12 18 19])&ismember("session_type",{''feedbackdelay'',''gonogo''})']);  % all NB cells
nonChAT = setdiff(nonChAT,{'n029_120222a_7.1' 'n045_130102a_4.2' 'n046_121229a_1.2' 'n046_130102a_6.3'});
% the first of the above cells excluded becuase of contamination from a
% light-activated unit; the other 3 because they have <=2 spikes
Hindex1 = getvalue('Hindex',nonChAT);
Hindex2 = getvalue('Hindex',ChAT);
Hindex2(13) = getvalue('Hindex','n046_130108b_4.2');   % use alternative session with better tagging session
Hindex = [Hindex1; Hindex2];
figure;
hist(Hindex,150)
fnm = [resdir 'Hindex3.fig'];
saveas(gcf,fnm)

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running TAGGINGSUMMARY...');  % progress indicator
global WB
WB(end+1) = wb;

% Reliability, latency, jitter
[Reliability Latency Jitter] = deal(nan(1,NumChAT));
SpikeNumberDistribution = cell(1,NumChAT);
cellids = ChAT;
for iC = 1:NumChAT
    cellid = cellids{iC};
    
    % Reliability, latency, jitter
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','stim','event','PulseOn','window',[-0.5 0.1],'isadaptive',2,...
        'baselinewin',[-0.005 0],'testwin',[0 0.01],'relative_threshold',0.03,...
        'jitterdefinition','burst','display',true,'rasterwindow',[-0.005 0.01]);
    figure(H.H_raster)
    xlim([-0.005 0.01])
    figure(H.H_psth)
    xlim([-0.005 0.01])
    Reliability(iC) = reliability;
    Latency(iC) = latency;
    Jitter(iC) = jitter;
    SpikeNumberDistribution{iC} = spikenumberdistribution;
    
    % Save
    if issave
        cellidt = regexprep(cellid,'\.','_');
        save([resdir cellidt '_TAGSUM.mat'],...
            'reliability','latency','jitter','spikenumberdistribution',...
            'B','M','lim1','lim2')
        fnm = [resdir cellidt '_TAGSUM_PSTH.jpg'];
        saveas(H.H_psth,fnm)
        fnm = [resdir cellidt '_TAGSUM_PSTH.fig'];
        saveas(H.H_psth,fnm)
        fnm = [resdir cellidt '_TAGSUM_RASTER.jpg'];
        set(H.H_raster,'PaperPositionMode','auto')
        set(H.H_raster,'InvertHardcopy','off')
        print(H.H_raster,'-djpeg',fnm)
        fnm = [resdir cellidt '_TAGSUM_RASTER.fig'];
        saveas(H.H_raster,fnm)
    end
    close(H.H_psth)
    close(H.H_raster)
    waitbar(iC/NumChAT)
end
close(wb)
if issave
    save([resdir 'tagging_latency_jitter.mat'],...
        'Reliability','Latency','Jitter','SpikeNumberDistribution',...
        'lim1','lim2','cellids')
end

keyboard

% Latency CDF
blue = [0 153 255] / 255;
figure
edges = 0:0.0001:max(Latency)+0.01;   % bin edges
dist_latency = histc(Latency,edges);   % latency histogram
dist_latency = [0 dist_latency(1:end-1)];   % values corresponding to the edges
dist_latency = dist_latency / sum(dist_latency);   % normalize
stairs(edges,cumsum(dist_latency),'Color',blue,'LineWidth',2)
axis tight
set(gca,'box','off','FontSize',12,'TickDir','out')
xlabel('Latency')
setmyplot_balazs

% Jitter CDF
figure
edges = 0:0.0001:max(Jitter)+0.01;   % bin edges
dist_jitter = histc(Jitter,edges);   % jitter histogram
dist_jitter = [0 dist_jitter(1:end-1)];   % values corresponding to the edges
dist_jitter = dist_jitter / sum(dist_jitter);   % normalize
stairs(edges,cumsum(dist_jitter),'Color',blue,'LineWidth',2)
axis tight
set(gca,'box','off','FontSize',12,'TickDir','out')
xlabel('Jitter')
setmyplot_balazs