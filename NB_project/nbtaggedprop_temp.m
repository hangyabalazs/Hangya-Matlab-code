function nbtaggedprop
%NBTAGGEDPROP   Properties of putative tagged neurons.

% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Find putative tagged cells
Lratio = getvalue('Lr_PC');
ID = getvalue('ID_PC');
Hindex = getvalue('Hindex');
R = getvalue('R');
ptinx = ID > 20 & Lratio < 0.15 & Hindex < 0.01 & R > 0.9;
putative_tagged = CELLIDLIST(ptinx);

% H-index, reliability, latency, jitter conditioned on burst type
NumCells = length(putative_tagged)';
for k = 1:NumCells
    cellid = putative_tagged{k};
    
    SE = loadcb(cellid,'StimEvents');
    burst_types = sort(unique(SE.BurstNPulse),'ascend');
    burst_types(isnan(burst_types)) = [];
    NumBurstTypes = length(burst_types);
    
    % Plot light-evoked and spont. spikes
    HS = plotwaveforms(cellid,'correlation',true,'maxnum',30;
    A1 = HS.activation_start;   % start of the detected stimulation period
    A2 = HS.activation_end;   % end of the detected stimulation period
    
    % Efficiency, latency and jitter for 'PulseOn'
    [E_burston L_burston J_burston B_burston M_burston] = ...
        reliability_latency_jitter(cellid,'event','PulseOn',...
        'activation_start',A1,'activation_end',Aend);
    
    % Tagging index for 'BurstOn'
    [Hindex_burston D_KL_burston] = tagging_index(cellid,'event','BurstOn');
    
    % Efficiency, latency and jitter for 'BurstOn'
    [E_burston L_burston J_burston B_burston M_burston] = ...
        reliability_latency_jitter(cellid,'event','BurstOn',...
        'activation_start',A1,'activation_end',Aend);
    
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
        [Hindex_frequency(bt) D_KL_frequency(bt) B_frequency(bt) M_frequency(bt)] = ...
            tagging_index(cellid,'event_filter','BurstNPulse_maxPower','filterinput',fi);  % H-index, D_KL
        [E_frequency(bt) L_frequency(bt) J_frequency(bt)] = reliability_latency_jitter(cellid,...
            'event_filter','BurstNPulse_maxPower','filterinput',fi,...
            'activation_start',A1,'activation_end',Aend);  % efficiency, latency, jitter
        
    end

    % Plot MClust projections
    plot_mclust_projections(cellid)
end




% Distance from light-evoked noise


% BurstOn and PulseOn PSTH