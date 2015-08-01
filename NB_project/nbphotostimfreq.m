function nbphotostimfreq
%NBPHOTOSTIMFFREQ   Light-response reliability.
%   NBPHOTOSTIMFFREQ calculates the reliability of light-response as a
%   function of stimulation frequency. Evoked spikes are auto-detected
%   based on doubly adaptive light-evoked PSTHs (see TAGGINGSUMMARY).
%
%   See also TAGGINGSUMMARY.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   3-Feb-2014

%   Edit log: BH 2/3/14

% Directories
global DATAPATH
resdir = [DATAPATH 'NB\taggingsummary\reliability2\'];
issave = true;

% Cells
selstr = ['"ChAT+"~=0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
cellids = ChAT;
cellids(11) = [];   % n046_121229x_1.1 - half session w/o stimulation
% 'n046_121230a_1.2' includes freely moving tagging session only to
% avoid drift effect! (valid_trials manually edited)
% 'n045_121217x_4.6' includes only after-behav. tagging session (cell could
% only be clustered from the second half of rec. - valid_trials manually
% edited)
NumCells = length(cellids);

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running NBPHOTOSTIMFREQ...');  % progress indicator
global WB
WB(end+1) = wb;

% Reliability
for iC = 1:NumCells
    main(cellids{iC},resdir,issave)
    waitbar(iC/NumCells)
end
close(wb)

% -------------------------------------------------------------------------
function main(cellid,resdir,issave)

% Light bursts
SE = loadcb(cellid,'StimEvents');
burst_types = sort(unique(SE.BurstNPulse),'ascend');
burst_types(isnan(burst_types)) = [];
NumBurstTypes = length(burst_types);

% Calculate the same tagging variables for bursts of different frequencies
[Reliability Latency Jitter] = deal(nan(1,NumBurstTypes));
SpikeNumberDistribution = cell(1,NumBurstTypes);
cellidt = regexprep(cellid,'\.','_');
for bt = 1:NumBurstTypes
    
    % Reliability, latency, jitter
    fi = struct('BurstNPulse',burst_types(bt));
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','stim','event','PulseOn','event_filter','BurstNPulse_maxPower',...
        'filterinput',fi,'window',[-0.5 0.1],'isadaptive',2,...
        'baselinewin',[-0.005 0],'testwin',[0 0.01],'relative_threshold',0.03,...
        'jitterdefinition','burst','display',true,'rasterwindow',[-0.005 0.01]);
    figure(H.H_raster)
    xlim([-0.005 0.01])
    figure(H.H_psth)
    xlim([-0.005 0.01])
    Reliability(bt) = reliability;
    Latency(bt) = latency;
    Jitter(bt) = jitter;
    SpikeNumberDistribution{bt} = spikenumberdistribution;
    
    % Save
    if issave
        fnm = [resdir cellidt '_' num2str(bt) '_TAGSUM_PSTH.jpg'];
        saveas(H.H_psth,fnm)
        fnm = [resdir cellidt '_' num2str(bt) '_TAGSUM_PSTH.fig'];
        saveas(H.H_psth,fnm)
        fnm = [resdir cellidt '_' num2str(bt) '_TAGSUM_RASTER.jpg'];
        set(H.H_raster,'PaperPositionMode','auto')
        set(H.H_raster,'InvertHardcopy','off')
        print(H.H_raster,'-djpeg',fnm)
        fnm = [resdir cellidt '_' num2str(bt) '_TAGSUM_RASTER.fig'];
        saveas(H.H_raster,fnm)
    end
    close(H.H_psth)
    close(H.H_raster)
end
if issave
    save([resdir cellidt '_tagging_reliability.mat'],...
        'Reliability','Latency','Jitter','SpikeNumberDistribution',...
        'lim1','lim2','burst_types')
end