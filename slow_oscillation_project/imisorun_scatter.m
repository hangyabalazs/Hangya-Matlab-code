function imisorun_scatter
%IMISORUN_SCATTER   Correlations among slow oscillation propagation properties.
%   IMISORUN_SCATTER calculates linear correlation among the following
%   parameters: propagation time, speed, distance and connection
%   strengrth.
%
%   See also IMISORUN.

% Directories
global DATADIR
global DATAPATH
pat = 'oiti40_mojzsis';
eg = '46';
inpdir = [DATAPATH 'Ulbert\OITI_40_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_40_EEG_' eg '\Scatter'];
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
nm_cols = 5;           % number of columns on the electrode grid (4 or 5)
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

% Scatter plots
H = figure;      % propagation time vs. connection strength
rIMn = rIM(rIM>0);
rIMLn = rIML(rIML>0);
plot(rIMn,rIMLn,'.')
xlabel('connection strength (Imax)')
ylabel('propagation time (max. loc.)')
[b,bint,r,rint,stats] = regress(rIMn,[ones(length(rIMLn),1) rIMLn]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(rIMn,rIMLn);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
fns = ['time_strength_' eg '.fig'];
saveas(H,fns)

H = figure;      % propagation speed vs. connection strength 
dmm = dm / 100;     % in meters
rIMLs = rIML / sr;  % in seconds
V = repmat(dmm,[1 1 tmax]) ./ rIMLs;
V = V(~isnan(V)&~isinf(V)&(V>0));
plot(rIMn,V,'.')
xlabel('connection strength (Imax)')
ylabel('propagation speed')
[b,bint,r,rint,stats] = regress(rIMn,[ones(length(V),1) V]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(rIMn,V);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
fns = ['speed_strength_' eg '.fig'];
saveas(H,fns)

H = figure;      % propagation time vs. propagation distance
ar = Adj .* repmat(dm,[1 1 tmax]);   % adjacency matrix with distances
arn = ar(ar>0);
plot(arn,rIMLn,'.')
xlabel('propagation distance (arrow length)')
ylabel('propagation time (max. loc.)')
[b,bint,r,rint,stats] = regress(arn,[ones(length(rIMLn),1) rIMLn]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(arn,rIMLn);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
fns = ['time_distance_' eg '.fig'];
saveas(H,fns)

H = figure;      % propagation speed vs. propagation distance
plot(arn,V,'.')
xlabel('propagation distance (arrow length)')
ylabel('propagation speed')
[b,bint,r,rint,stats] = regress(arn,[ones(length(V),1) V]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(arn,V);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
fns = ['speed_distance_' eg '.fig'];
saveas(H,fns)

H = figure;      % connection strength vs. propagation distance
plot(arn,rIMn,'.')
xlabel('propagation distance (arrow length)')
ylabel('connection strength (Imax)')
[b,bint,r,rint,stats] = regress(arn,[ones(length(rIMn),1) rIMn]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(arn,rIMn);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
fns = ['strength_distance_' eg '.fig'];
saveas(H,fns)



% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end