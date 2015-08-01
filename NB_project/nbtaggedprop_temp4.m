function error_list = nbtaggedprop(I,issave)
%NBTAGGEDPROP   Properties of putative tagged neurons.
%   NBTAGGEDPROP(I,ISSAVE) performs analysis aiding the definite decision
%   about tagging. Input parameters: 
%       I - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, putative tagged
%           cells are selected (ID>20, L-ratio<0.15, H index<0.01, R>0.9;
%           see NBTAGGING, LRATIO, NBISSTIM and NBSPIKESHAPECORR for
%           details on these measures)
%       ISSAVE - controls saving
%
%   The following analyses are performed (output variables saved in mat
%   files and figures in pdf):
%       Reliability, latency and jitter of evoked spikes after 'BurstOn',
%           'PulseOn' and frequency-restricted 'PulseOn' events; see
%           RELIABILITY_LATENCY_JITTER for details
%       H index for 'BurstOn', 'PulseOn' and frequency-restricted 'PulseOn'
%           events; see TAGGING_INDEX for details
%       Waveforms of spontaneous and light-evoked spikes; see PLOTWAVEFORMS
%           for details
%       All projections of feature data in the Energy-Amplitue space with
%           the putative tagged cells shown in orange and the light-evoked 
%           spikes overlayed in blue and ; see PLOT_MCLUST_PROJECTIONS2 for
%           details
%       Cluster quality indices restricted to light-evoked spikes in the
%           Energy-Amplitude as well as in the Energy-WavePC1 space; see
%           LIGHTCLUSTER and LRATIO for details
%       Raster plot and PSTH aligned to 'BurstOn' and 'PulseOn' events
%           (only 1000 pulses shown in the latter); see PLOT_RASTER_PSTH,
%           VIEWCELL2B and PLOT_RASTER2A for details.
%
%   ERROR_LIST = NBTAGGEDPROP(I,ISSAVE) returns a structure with all caught
%   errors. ERROR_LIST includes cell IDs, error messages and the captured
%   exception variables. A time stamped ERROR_LIST is also assigned in base
%   workspace and saved to the results directory automatically.
%
%   See also NBTAGGING, LRATIO, NBISSTIM, NBSPIKESHAPECORR,
%   RELIABILITY_LATENCY_JITTER, TAGGING_INDEX, PLOTWAVEFORMS,
%   PLOT_MCLUST_PROJECTIONS2, LIGHTCLUSTER, PLOT_RASTER_PSTH,
%   VIEWCELL2B and PLOT_RASTER2A.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   05-Oct-2012

%   Edit log: BH 5/10/12, 4/19/13

% Pass the control to the user in case of error
dbstop if error
 
% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1
    I = [];
end
 
% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'taggedprop3_new' fs];
 
% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Find putative tagged cells
if isempty(I)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    Hindex = getvalue('Hindex');
    R = getvalue('R');
    ptinx = ID > 20 & Lratio < 0.15 & Hindex < 0.01 & R > 0.9;
    I = find(ptinx);
    I = I(I>2309);
    putative_tagged = CELLIDLIST(I);
else
    Hindex = getvalue('Hindex');
    if isnumeric(I)
        putative_tagged = CELLIDLIST(I);
    else
        putative_tagged = I;
    end
end

% All the analysis promised in the help
NumCells = length(putative_tagged);
error_list = struct('cellid',{},'message',{},'exception',{});  % keep a list of caught errors 
errinx = 0;
for k = 1:NumCells
    cellid = putative_tagged{k};
    try
        main(cellid,issave,resdir)  % every analysis
    catch ME
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
        errinx = errinx + 1;  % error counter
        error_list(errinx).cellid = cellid;   % record problematic cell ID
        error_list(errinx).message = ME.message;   % error message
        error_list(errinx).exception = ME;   % exception structure
    end
end

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
list_name = ['error_list_' dsr];  % variable name
list_name = regexprep(list_name,'-','_');
list_name = regexprep(list_name,' ','_');
assignin('base',list_name,error_list)   % assign error list in base workspace
error_fname = fullfile(resdir,[list_name '.mat']);   % file name
save(error_fname,'error_list')   % save error list

% -------------------------------------------------------------------------
function main(cellid,issave,resdir)

% Add 'PulseOn' event if missing
ST = loadcb(cellid,'STIMSPIKES');
if isequal(findcellstr(ST.events(:,1),'PulseOnS'),0)
    [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
    prealignSpikes(cellid,'events',...
        {'PulseOnS' 'PulseOn' 'PulseOn' [lim1 lim2]},...
        'epochs',[],'filetype','stim','ifsave',1,'ifappend',1)
else
    [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
    prealignSpikes(cellid,'events',...
        {'PulseOnS' 'PulseOn' 'PulseOn' [lim1 lim2]},...
        'epochs',[],'filetype','stim','ifsave',1,'writing_behavior','replace')
end

SE = loadcb(cellid,'StimEvents');
burst_types = sort(unique(SE.BurstNPulse),'ascend');
burst_types(isnan(burst_types)) = [];
NumBurstTypes = length(burst_types);

% Efficiency, latency and jitter for 'PulseOn'
[E_pulseon L_pulseon J_pulseon B_pulseon M_pulseon A1 A2] = ...
    reliability_latency_jitter(cellid,'event','PulseOn');

% H-index for 'PulseOn'
Hindex_pulseon = HPO;

% Efficiency, latency and jitter for 'BurstOn'
[E_burston L_burston J_burston B_burston M_burston] = ...
    reliability_latency_jitter(cellid,'event','BurstOn');

% H-index for 'BurstOn'
[Hindex_burston D_KL_burston] = tagging_index(cellid,'event','BurstOn');  %#ok<NASGU> % H-index, D_KL

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
HE = figure('Position',[624 110 900 868]);
S = set_subplots(4,1,0.065,0.065);
axes(S(1)) %#ok<*MAXES>
bar(S(1),E_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
set(S(1),'XTickLabel',num2cell(burst_types/2))
L1 = line(xlim,[E_pulseon E_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
L2 = line(xlim,[E_burston E_burston],'Color',[255 204 0]/255,'LineWidth',2);
legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
title('Efficiency')

% Plot latency and jitter
axes(S(2))
bar(S(2),L_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
hold on
errorbar(L_frequency,J_frequency,'k+')
set(S(2),'XTickLabel',num2cell(burst_types/2))
L1 = line(xlim,[L_pulseon L_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
line(xlim,[L_pulseon+J_pulseon L_pulseon+J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
line(xlim,[L_pulseon-J_pulseon L_pulseon-J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
L2 = line(xlim,[L_burston L_burston],'Color',[255 204 0]/255,'LineWidth',2);
line(xlim,[L_burston+J_burston L_burston+J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
line(xlim,[L_burston-J_burston L_burston-J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
title('Latency')

% Plot evoked firing rate
axes(S(3))
bar(S(3),M_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
set(S(3),'XTickLabel',num2cell(burst_types/2))
L1 = line(xlim,[M_pulseon M_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
L2 = line(xlim,[M_burston M_burston],'Color',[255 204 0]/255,'LineWidth',2);
L3 = line(xlim,[B_pulseon B_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
L4 = line(xlim,[B_burston B_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
legend([L1 L2 L3 L4],{'PulseOn' 'BurstOn' 'Baseline' 'Baseline'},...
    'Location','EastOutside')
title('Evoked firing rate')

% Plot H-index
axes(S(4))
bar(S(4),Hindex_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
set(S(4),'XTickLabel',num2cell(burst_types/2))
L1 = line(xlim,[Hindex_pulseon Hindex_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
L2 = line(xlim,[Hindex_burston Hindex_burston],'Color',[255 204 0]/255,'LineWidth',2);
line(xlim,[0.01 0.01],'Color','r','LineWidth',2)
legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
title('H index')

% Plot light-evoked and spont. spikes
HW = plotwaveforms(cellid,'correlation',true,'maxnum',5000);

% Plot MClust projections
HM = plot_mclust_projections2(cellid,'feature_names',{'Amplitude','Energy'},...
    'stim_period',[A1 A2]);

% Distance from light-evoked noise
[lID_amp lLr_amp Pref Pref2] = lightcluster(cellid,'feature_names',{'Amplitude','Energy'},...
    'stim_period',[A1 A2]);
[lID_PC lLr_PC] = lightcluster(cellid,'feature_names',{'Amplitude','WavePC1'},...
    'stim_period',[A1 A2]);
HD = figure;
axes;
str = {'Cluster quality measures for light-evoked spikes:';...
    ' ';...
    ['ID (Amplitude, Energy): ' num2str(lID_amp)];...
    ['L-ratio (Amplitude, Energy): ' num2str(lLr_amp)];...
    ['ID (WavePC1, Energy): ' num2str(lID_PC)];...
    ['L-ratio (WavePC1, Energy): ' num2str(lLr_PC)];...
    ['preference: ' num2str(Pref)];...
    ['preference: ' num2str(Pref2)]};
uicontrol('Style','text','Unit','normalized','Position',...
    [0.18 0.3 0.65 0.5],'FontSize',12,'HorizontalAlignment',...
    'left','String',str,'BackgroundColor',[1 1 1]);
axis off

% BurstOn and PulseOn PSTH
HR = plot_raster_psth(cellid,'BurstOn',true,'PulseOn',true);

% Save
if issave
    save([resdir 'TAGGEDPROP_' regexprep(cellid,'\.','_') '.mat'],...
        'Hindex_burston','D_KL_burston',...
        'E_pulseon','L_pulseon','J_pulseon','B_pulseon','M_pulseon','A1','A2',...
        'E_burston','L_burston','J_burston','B_burston','M_burston',...
        'E_frequency','L_frequency','J_frequency','B_frequency','M_frequency',...
        'lLr_amp','lID_amp','lLr_PC','lID_PC','Pref','Pref2')
    
    % Write figures to pdf
    pdfname = fullfile(resdir,['TAGGEDPROP_' regexprep(cellid,'\.','_') '.pdf']);
    writefigs(HR,pdfname)
    writefigs(HE,pdfname)
    writefigs(HW,pdfname)
    writefigs(HM,pdfname)
    writefigs(HD,pdfname)
    close all
end

% -------------------------------------------------------------------------
function writefigs(H,pdfname)
% Write figures to pdf
 
% Append to pdf
if isstruct(H)  % H is either a struct with figure handles, or a single fig. handle
    fls = fieldnames(H);
    for fs = 1:length(fls)
        h = H.(fls{fs});
        if ishandle(h)
            export_fig(h,'-append',pdfname, '-zbuffer');  % write to pdf
        end
    end
else
    export_fig(H,'-append',pdfname, '-zbuffer');  % write to pdf
end