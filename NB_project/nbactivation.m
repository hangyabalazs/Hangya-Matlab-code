function nbactivation
%NBACTIVATION   Quantify spike timing.
%   NBACTIVATION calculates latency and jitter of spikes evoked by
%   punishment for cholinergic neurons. NBACTIVATION calls
%   RELIBILITY_LATENCY_JITTER, which in turn calls ULTIMATE_PSTH (see
%   further details therein). Only the first spikes in the trials are used
%   for jitter calculation, assuming burst firing. The resulting variables
%   are saved along with a PSTH and a raster plot.
%
%   See also RELIBILITY_LATENCY_JITTER and ULTIMATE_PSTH.

% Directories
global DATAPATH
responsetype = 'Hit';
resdir = [DATAPATH 'NB\' responsetype '_activation_newdata\'];
issave = true;

% Cholinergic neurons
ChAT1 = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
ChAT1 = [ChAT1 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
pChAT1 = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT1 pChAT1];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running NBACTIVATION...');  % progress indicator
global WB
WB(end+1) = wb;

% Reliability, latency, jitter
[Reliability Latency Jitter] = deal(nan(1,NumChAT));
SpikeNumberDistribution = cell(1,NumChAT);
cellids = allChAT;
for iC = 1:NumChAT
    cellid = cellids{iC};
    
    % Checking whether 'DeliverFeedback' event is available
    sesstype = getvalue('session_type',cellid);
    if isequal(sesstype,{'feedbackdelay'})
        alignevent = 'DeliverFeedback';
    else
        switch responsetype
            case 'FA'
                alignevent = 'LeftPortIn';
            case 'Hit'
                alignevent = 'LeftWaterValveOn';
        end
    end
    
    % Filter for trials
    switch responsetype
        case 'FA'
            trialfilter = 'FalseAlarm==1';
        case 'Hit'
            trialfilter = 'Hit==1';
    end
    
    % Reliability, latency, jitter
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
        'event_filter','custom','filterinput',trialfilter,'isadaptive',2,...
        'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
        'jitterdefinition','burst','display',true);
    Reliability(iC) = reliability;
    Latency(iC) = latency;
    Jitter(iC) = jitter;
    SpikeNumberDistribution{iC} = spikenumberdistribution;
    
    % Save
    if issave
        cellidt = regexprep(cellid,'\.','_');
        save([resdir cellidt '_' responsetype '_ACTIVATION.mat'],...
            'reliability','latency','jitter','spikenumberdistribution',...
            'B','M','lim1','lim2')
        fnm = [resdir cellidt '_' responsetype '_ACTIVATION_PSTH.jpg'];
        saveas(H.H_psth,fnm)
        fnm = [resdir cellidt '_' responsetype '_ACTIVATION_PSTH.fig'];
        saveas(H.H_psth,fnm)
        fnm = [resdir cellidt '_' responsetype '_ACTIVATION_RASTER.jpg'];
        set(H.H_raster,'PaperPositionMode','auto')
        set(H.H_raster,'InvertHardcopy','off')
        print(H.H_raster,'-djpeg',fnm)
        fnm = [resdir cellidt '_' responsetype '_ACTIVATION_RASTER.fig'];
        saveas(H.H_raster,fnm)
    end
    close(H.H_psth)
    close(H.H_raster)
    waitbar(iC/NumChAT)
end
close(wb)
save([resdir responsetype '_activation_timing.mat'],...
    'Reliability','Latency','Jitter','SpikeNumberDistribution','cellids')