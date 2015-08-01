function iwhatpropagates
%IWHATPROPAGATES   Slow oscillation propagation.
%   IWHATPROPAGATES averages data windows where significant MI-shift
%   occurred. Data is plotted and saved for each channels.
%
%   See also IMISORUN.

% Directories
global DATADIR
global DATAPATH
pat = 'oitin1_wittner';
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
eg = '157';
esg = [1900 1960];
inpdir = [DATAPATH 'Ulbert\OITI_N1_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_N1_EEG_' eg '\WhatPropagates\'];
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end

% Load MI map
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)

% Significance level
fn = [inpdir 'siglev_EEG' eg];
load(fn)
siglev = sl(4,2);   % sig. lev.: 0.0001

% Transform input variables
Adj = rIMax;        % adjacency matrix
Adj(rIMax<siglev) = NaN;
Adj(isnan(Adj)) = 0;
Adj(Adj>0) = 1;     % time-varying adjacency matrix
rIM = rIMax;
rIM(rIMax<siglev) = NaN;
rIM(isnan(rIM)) = 0;
rIML = rIMaxLoc;
rIML(rIMax<siglev) = NaN;
rIML(isnan(rIML)) = 0;

% Recording properties
chnum = size(Adj,1);   % number of grid electrodes
tmax = size(Adj,3);    % time axis length
nm_rows = 4;           % number of rows on the electrode grid
nm_cols = 4;           % number of columns on the electrode grid (4 or 5)
if ~isequal(nm_rows*nm_cols,chnum)
    error('Improper number for grid electrodes.')
end
sr = 1000;             % sampling rate

% Load raw data
ddir = [DATADIR '\human_SO\' pat '\grid\mat_EEG_' eg '\'];
fnm = ['EEG_' eg '_' num2str(esg(1)) '_' num2str(esg(2)) '_rs_filt01_40.mat'];
load([ddir fnm])

% What propagates?
stck = cell(1,chnum);
for i = 1:chnum
    for k = 1:tmax
        if any(rIM(i,:,k))
            stck{i}(1:1000,end+1) = data((k-1)*100+1:(k-1)*100+1000,i);
        end
    end
end

% Plot
H = figure;
for i = 1:chnum
    subplot(nm_rows,nm_cols,i);
    hold on
    if ~isempty(stck{i})
        plot(mean(stck{i},2))
    end
end

% Save
saveas(H,['whatpropagates_' eg])