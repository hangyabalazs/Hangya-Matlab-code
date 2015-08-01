%VIPMAKEPDF   Writes pdf files.
%   VIPMAKEPDF writes pdf files for all cell significantly influenced by
%   light stimulation. Light-pulse aligned raster plot, SDF as well as
%   auto-correlation and spike shape is saved to pdf.

%% load

load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')  % load all variables fo light-effects

%% groups

% Group cells based on light effects
tagged = logical(vldty) & (isact==2);   % tagged group
inx_act = logical(vldty) & (isact==1);   % indices for activated cells
inx_inh = logical(vldty) & (isinh) & (baseline>1);   % indices for inhibited cells - only if baseline rate > 1 Hz
activated = setdiff(find(inx_act&~inx_inh),tagged);  % activated and not inhibited cells
inhibited = find(inx_inh&~inx_act);  % inhibited and not activated cells
ai = find(inx_act&inx_inh);  % cells that are both activated and inhibited
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));  % first activated then inhibited
inhibited_activated = ai(~inx);  % first inhibited then activated
activated = [activated; activated_inhibited];  % first activated, either inhibited or not
inhibited = [inhibited; inhibited_activated];  % first inhibited, either activated or not
tagged = find(tagged);  % indices for tagged cells

%% load figures and write pdf

% Data directories
global DATAPATH
inpdir1 = [DATAPATH 'VIP\vipisinfluenced11_A1_newcells3_p\'];   % for PSTHs and raster plots
inpdir2 = [DATAPATH 'VIP\autocorr_ACx\'];   % for auto-correlations
inpdir3 = [DATAPATH 'VIP\waveform_ACx\'];   % for spike shape
resdir1 = [DATAPATH 'VIP\A1_lighteffects_activated\'];
resdir2 = [DATAPATH 'VIP\A1_lighteffects_inhibited\'];

% Load CellBase
loadcb

% Activated cells - load figures
for k = 1:length(activated)
    cellid = CELLIDLIST{activated(k)};   % current CellID
    
    fn = [inpdir2 num2str(activated(k)) '_' cellid '_acorr.fig'];   % auto-correlation
    Ha = open(fn);
    
    fn = [inpdir3 num2str(activated(k)) '_' cellid '_wv.fig'];   % spike shape
    Hw = open(fn);
    
    fn = [inpdir1 'SDF_' cellid '.fig'];   % PSTH (SDF)
    Hs = open(fn);
    
    fn = [inpdir1 'RASTER_' cellid '.fig'];   % raster plot
    Hr = open(fn);
    
    pdfname = [resdir1 'LIGHTEFFECT_' regexprep(cellid,'\.','_') '.pdf'];
    writefigs(Hs,pdfname)
    writefigs(Hr,pdfname)
    writefigs(Hw,pdfname)
    writefigs(Ha,pdfname)
    close all
    
end

% Inhibited cells - load figures
for k = 1:length(inhibited)
    cellid = CELLIDLIST{inhibited(k)};   % current CellID
    
    fn = [inpdir2 num2str(inhibited(k)) '_' cellid '_acorr.fig'];   % auto-correlation
    Ha = open(fn);
    
    fn = [inpdir3 num2str(inhibited(k)) '_' cellid '_wv.fig'];   % spike shape
    Hw = open(fn);
    
    fn = [inpdir1 'SDF_' cellid '.fig'];   % PSTH (SDF)
    Hs = open(fn);
    
    fn = [inpdir1 'RASTER_' cellid '.fig'];   % raster plot
    Hr = open(fn);
    
    pdfname = [resdir2 'LIGHTEFFECT_' regexprep(cellid,'\.','_') '.pdf'];
    writefigs(Hs,pdfname)
    writefigs(Hr,pdfname)
    writefigs(Hw,pdfname)
    writefigs(Ha,pdfname)
    close all
    
end

%% load figures and write pdf for individual cells

% Data directories
global DATAPATH
inpdir1 = [DATAPATH 'VIP\vipisinfluenced11_A1_newcells3_p\'];   % for PSTHs and raster plots
inpdir2 = [DATAPATH 'VIP\autocorr_ACx\'];   % for auto-correlations
inpdir3 = [DATAPATH 'VIP\waveform_ACx\'];   % for spike shape
resdir1 = [DATAPATH 'VIP\A1_lighteffects_plus\'];

% Load CellBase
loadcb

% Load figures
k = 302;
cellid = CELLIDLIST{k};   % current CellID

fn = [inpdir2 num2str(k) '_' cellid '_acorr.fig'];   % auto-correlation
Ha = open(fn);

fn = [inpdir3 num2str(k) '_' cellid '_wv.fig'];   % spike shape
Hw = open(fn);

fn = [inpdir1 'SDF_' cellid '.fig'];   % PSTH (SDF)
Hs = open(fn);
    
    fn = [inpdir1 'RASTER_' cellid '.fig'];   % raster plot
    Hr = open(fn);
    
    pdfname = [resdir1 'LIGHTEFFECT_' regexprep(cellid,'\.','_') '.pdf'];
    writefigs(Hs,pdfname)
    writefigs(Hr,pdfname)
    writefigs(Hw,pdfname)
    writefigs(Ha,pdfname)
    close all