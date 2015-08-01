function nbtaggedprop(I,issave)
%NBTAGGEDPROP   Properties of putative tagged neurons.

% Pass the control to the user in case of error
dbstop if error

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = false;
end
if nargin < 1
    I = [];
end

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'taggedprop' fs];

% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Find putative tagged cells
if isempty(I)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    Hindex = getvalue('Hindex');
    R = getvalue('R');
    ptinx = ID > 20 & Lratio < 0.15 & Hindex < 0.01 & R > 0.9;
    fptinx = find(ptinx);
    putative_tagged = CELLIDLIST(ptinx);
else
    putative_tagged = CELLIDLIST(I);
end

% H-index, reliability, latency, jitter conditioned on burst type
NumCells = length(putative_tagged);
for k = 1:NumCells
    cellid = putative_tagged{k};
    
    SE = loadcb(cellid,'StimEvents');
    burst_types = sort(unique(SE.BurstNPulse),'ascend');
    burst_types(isnan(burst_types)) = [];
    NumBurstTypes = length(burst_types);
    
    % Efficiency, latency and jitter for 'PulseOn'
    [E_pulseon L_pulseon J_pulseon B_pulseon M_pulseon A1 A2] = ...
        reliability_latency_jitter(cellid,'event','PulseOn');
    
    % H-index for 'PulseOn'
    Hindex_pulseon = Hindex(fptinx(k));
    
    % Efficiency, latency and jitter for 'BurstOn'
    [E_burston L_burston J_burston B_burston M_burston] = ...
        reliability_latency_jitter(cellid,'event','BurstOn');
    
    % H-index for 'BurstOn'
    [Hindex_burston D_KL_burston] = tagging_index(cellid,'event','BurstOn');  %#ok<NASGU,ASGLU> % H-index, D_KL
    
    % Calculate the same tagging variables for bursts of different frequencies
    Hindex_frequency = nan(1,NumBurstTypes);
    D_KL_frequency = nan(1,NumBurstTypes);
    E_frequency = nan(1,NumBurstTypes);
    L_frequency = nan(1,NumBurstTypes);
    J_frequency = nan(1,NumBurstTypes);
    B_frequency = nan(1,NumBurstTypes);
    M_frequency = nan(1,NumBurstTypes);
    for bt = 1:NumBurstTypes
        
        fi = struct('BurstNPulse',burst_types(bt));
        [Hindex_frequency(bt) D_KL_frequency(bt)] = ...
            tagging_index(cellid,'event_filter','BurstNPulse_maxPower','filterinput',fi);  % H-index, D_KL
        [E_frequency(bt) L_frequency(bt) J_frequency(bt) B_frequency(bt) M_frequency(bt)] = reliability_latency_jitter(cellid,...
            'event_filter','BurstNPulse_maxPower','filterinput',fi);  % efficiency, latency, jitter
        
    end
    
    % Plot efficiency
    HE = figure('Position',[624 110 672 868]);
    S = set_subplots(4,1,0.05,0.05);
    axes(S(1)) %#ok<LAXES>
    bar(S(1),E_frequency,0.5,'EdgeColor','k','FaceColor','w')
    set(S(1),'XTickLabel',num2cell(burst_types/2))
    L1 = line(xlim,[E_pulseon E_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
    L2 = line(xlim,[E_burston E_burston],'Color',[255 204 0]/255,'LineWidth',2);
    legend([L1 L2],{'PulseOn' 'BurstOn'})
    title('Efficiency')
    
    % Plot latency and jitter
    axes(S(2)) %#ok<LAXES>
    bar(S(2),L_frequency,0.5,'EdgeColor','k','FaceColor','w')
    errorbar(S(2),J_frequency,'k+')
    set(S(2),'XTickLabel',num2cell(burst_types/2))
    L1 = line(xlim,[L_pulseon L_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
    line(xlim,[L_pulseon+J_pulseon L_pulseon+J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle',':');
    line(xlim,[L_pulseon-J_pulseon L_pulseon-J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle',':');
    L2 = line(xlim,[L_burston L_burston],'Color',[255 204 0]/255,'LineWidth',2);
    line(xlim,[L_burston+J_burston L_burston+J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle',':');
    line(xlim,[L_burston-J_burston L_burston-J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle',':');
    legend([L1 L2],{'PulseOn' 'BurstOn'})
    title('Latency')
    
    % Plot evoked firing rate
    axes(S(3)) %#ok<LAXES>
    bar(S(3),M_frequency,0.5,'EdgeColor','k','FaceColor','w')
    set(S(3),'XTickLabel',num2cell(burst_types/2))
    L1 = line(xlim,[M_pulseon M_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
    L2 = line(xlim,[M_burston M_burston],'Color',[255 204 0]/255,'LineWidth',2);
    L3 = line(xlim,[B_pulseon B_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle',':');
    L4 = line(xlim,[B_burston B_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle',':');
    legend([L1 L2 L3 L4],{'PulseOn' 'BurstOn' 'Baseline (PulseOn)' 'Baseline (BurstOn)'})
    title('Evoked firing rate')
    
    % Plot H-index
    axes(S(4)) %#ok<LAXES>
    bar(S(4),H_frequency,0.5,'EdgeColor','k','FaceColor','w')
    set(S(4),'XTickLabel',num2cell(burst_types/2))
    L1 = line(xlim,[Hindex_pulseon Hindex_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
    L2 = line(xlim,[H_burston H_burston],'Color',[255 204 0]/255,'LineWidth',2);
    line(xlim,[0.01 0.01],'Color','r','LineWidth',2)
    legend([L1 L2],{'PulseOn' 'BurstOn'})
    title('H index')
    
    % Plot light-evoked and spont. spikes
    HW = plotwaveforms(cellid,'correlation',true,'maxnum',30);

    % Plot MClust projections
    HM = plot_mclust_projections2(cellid,'feature_names',{'Amplitude','Energy'},...
        'stim_period',[A1 A2]);
    
    % Distance from light-evoked noise
    [lID_amp lLr_amp] = lightcluster(cellid,'feature_names',{'Amplitude','Energy'},...
        'stim_period',[A1 A2]);
    [lID_PC lLr_PC Pref] = lightcluster(cellid,'feature_names',{'Amplitude','WavePC1'},...
        'stim_period',[A1 A2]);
    HD = figure;
    axes;
    str = {['ID (Amplitude, Energy): ' num2str(lID_amp)];...
        ['L-ratio (Amplitude, Energy): ' num2str(lLr_amp)];...
        ['ID (WavePC1, Energy): ' num2str(lID_PC)];...
        ['L-ratio (WavePC1, Energy): ' num2str(lLr_PC)];...
        ['preference: ' num2str(Pref)]};
    text(0,0,str)
    
    % BurstOn and PulseOn PSTH
    HR = plot_raster_psth(cellid,'BurstOn',true,'PulseOn',true);
    
    % Save
    if issave
        save([resdir 'TAGGEDPROP_' regexprep(cellid,'\.','_') '.mat'],...
            'Hindex_burston','D_KL_burston',...
            'E_pulseon','L_pulseon','J_pulseon','B_pulseon','M_pulseon','A1','A2',...
            'E_burston','L_burston','J_burston','B_burston','M_burston',...
            'E_frequency','L_frequency','J_frequency','B_frequency','M_frequency',...
            'lLr_amp','lID_amp','lLr_PC','lID_PC','Pref')
        
        % Write figures to pdf
        pdfname = fullfile(resdir,['TAGGEDPROP_' regexprep(cellid,'\.','_') '.pdf']);
        writefigs(HR,pdfname)
        writefigs(HE,pdfname)
        writefigs(HW,pdfname)
        writefigs(HM,pdfname)
        writefigs(HD,pdfname)
    end
end

% -------------------------------------------------------------------------
function writefigs(H,pdfname)
% Write figures to pdf

% Append to pdf
if isstruct(H)  % H is either a struct with figure handles, or a single fig. handle
    fls = fieldnames(H);
    for fs = 1:length(fls)
        export_fig(H.(fls),'-append',pdfname, '-zbuffer');  % write to pdf
    end
else
    export_fig(H,'-append',pdfname, '-zbuffer');  % write to pdf
end