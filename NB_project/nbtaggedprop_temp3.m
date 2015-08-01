function error_list = nbtaggedprop_temp3(I,issave)
%NBTAGGEDPROP   Properties of putative tagged neurons.
%   NBTAGGEDPROP(I,ISSAVE) performs analysis aiding the definite decision
%   about tagging. Input parameters: 
%       I - index set to CELLIDLIST (see CellBase documentation); if empty
%           or not specified, putative tagged cells are selected (ID>20,
%           L-ratio<0.15, H index<0.01, R>0.9; see NBTAGGING, LRATIO,
%           NBISSTIM and NBSPIKESHAPECORR for details on these measures)
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
%   See also NBTAGGING, LRATIO, NBISSTIM, NBSPIKESHAPECORR,
%   RELIABILITY_LATENCY_JITTER, TAGGING_INDEX, PLOTWAVEFORMS,
%   PLOT_MCLUST_PROJECTIONS2, LIGHTCLUSTER, PLOT_RASTER_PSTH,
%   VIEWCELL2B and PLOT_RASTER2A.

%   Edit log: BH 5/10/12
 
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
resdir = [DATAPATH 'NB' fs 'taggedprop3_new_temp' fs];
 
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
%     I = I(I>2309);
    putative_tagged = CELLIDLIST(I);
else
    Hindex = getvalue('Hindex');
    putative_tagged = CELLIDLIST(I);
end

% All the analysis promised in the help
NumCells = length(putative_tagged);
error_list = struct('cellid',{},'message',{},'exception',{});
errinx = 0;
for k = 50:NumCells
    cellid = putative_tagged{k};
    disp([num2str(k) '   ' cellid])
    try
        HPO = Hindex(I(k)); % already calculated by NBTAGGING, so pass it on
        main(cellid,issave,resdir,HPO)  % every analysis
    catch ME
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
        errinx = errinx + 1;
        error_list(errinx).cellid = cellid;
        error_list(errinx).message = ME.message;
        error_list(errinx).exception = ME;
    end
end

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
list_name = ['error_list_' dsr];  % variable name
list_name = regexprep(list_name,'-','_');
list_name = regexprep(list_name,' ','_');
assignin('base',list_name,error_list)   % assign error list in base workspace

% -------------------------------------------------------------------------
function main(cellid,issave,resdir,HPO)

% Efficiency, latency and jitter for 'PulseOn'
[E_pulseon L_pulseon J_pulseon B_pulseon M_pulseon A1 A2] = ...
    reliability_latency_jitter(cellid,'event','PulseOn');

% Distance from light-evoked noise
[lID_amp lLr_amp Pref Pref2] = lightcluster_temp(cellid,'feature_names',{'Amplitude','Energy'},...
    'stim_period',[A1 A2]);
[lID_PC lLr_PC] = lightcluster_temp(cellid,'feature_names',{'Amplitude','WavePC1'},...
    'stim_period',[A1 A2]);