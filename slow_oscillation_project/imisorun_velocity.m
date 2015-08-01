function imisorun_velocity
%IMISORUN_VELOCITY   Slow oscillation propagation speed.
%   IMISORUN_VELOCITY calculates the speed of slow oscillation propagation.
%
%   See also IMISORUN.

% Directories
global DATAPATH
eg = '40';
inpdir = [DATAPATH 'Ulbert\OITI_39_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_39_EEG_' eg '\Speed'];
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end

% Load
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
nm_cols = 4;           % number of columns on the electrode grid
if ~isequal(nm_rows*nm_cols,chnum)
    error('Improper number for grid electrodes.')
end
sr = 1000;             % sampling rate

% Distance matrix
dm = zeros(chnum,chnum);
for x = 1:chnum
    for y = 1:chnum
        [i j] = gridind2sub(x,nm_cols);
        [k l] = gridind2sub(y,nm_cols);
        dm(x,y) = sqrt((i-k)^2 + (j-l)^2);
    end
end

% Propagation speed
dmm = dm / 100;     % in meters
rIMLs = rIML / sr;  % in seconds
V = repmat(dmm,[1 1 tmax]) ./ rIMLs;
V = V(~isnan(V)&~isinf(V)&(V>0));

lims = [1/32 0.0625 0.125 0.25 0.5 1 2 4 8 16];
N = length(V(:));
hs(1,1) = length(V(V<lims(1))) / N;
for k = 2:length(lims)
    hs(1,k) = length(V(V>lims(k-1)&V<lims(k))) / N;
end
hs(1,length(lims)+1) = length(V(V>lims(end))) / N;
speed = hs;
allspeed = V(:);

% Normalized arrow length distribution
ar = Adj .* repmat(dm,[1 1 tmax]);   % adjacency matrix with distances
arn = ar(ar>0);
[nm xout] = hist(arn,50);      % arrow length distribution
[nm2 xout2] = hist(dm(dm>0),50);
if ~isequal(xout,xout2)
    error('Technical error 522.')
end
nmm = nm ./ nm2;    % normalization by the abundance of a given distance
nmm2 = nan2zero(nmm);
arrowlength = nmm2 / sum(nmm2);
arrowlength_orig = nm / sum(nm);
alx = xout;

% Connection strength distribution
rIMn = rIM(rIM>0);
[na xa] = hist(rIMn);
arrowstrength = na;
asx = xa;

save('vars.mat','allspeed','speed','arrowlength','arrowlength_orig','alx','arrowstrength','asx')


% figure;bar(hs,'stacked')
% saveas(H,['velocitydist_' eg])



% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end