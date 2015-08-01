function iwithinsubjcorr2(pno,pt,egs)
%IWITHINSUBJCORR2   Within-subject correlations.
%   IWITHINSUBJCORR2 calculates within-subject correlation of propagation
%   speed, distance and strength. See IMISORUN for details on the above
%   variables. Distribution distances are measured by symmetric
%   KL-divergence.
%
%   IWITHINSUBJCORR2(PNO,PT,EGS) requires 3 string input arguments:
%       PNO: patient ID
%       PT: patient name
%       EGS: EEG file IDs
%
%   See also IMISORUN and IMISORUN_MTX.

% Directories
global DATADIR
global DATAPATH
% pno = num2str(39);
% % pno = 'n1';
% pt = 'gaal'
% egs = [{'10'} {'14'} {'40'}];
pat = ['oiti' pno '_' pt];
resdir = [DATAPATH 'Ulbert\Summary_temp\'];

% Load
speed_dist = cell(1,3);
strength_dist = cell(1,3);
length_dist = cell(1,3);
normlength_dist = cell(1,3);
for sgs = 1:3
    eg = egs{sgs};
    inpdir = [DATAPATH 'Ulbert\OITI_' pno '_EEG_' eg '\Mtx\'];
    
    fn = [inpdir 'mtxs_' eg '.mat'];    % MI map
    load(fn)
    
    V = Speed;      % speed distribution
    lims = [1/32 0.0625 0.125 0.25 0.5 1 2 4 8 16];
    N = length(V(:));
    hs(1,1) = length(V(V<lims(1))) / N;
    for k = 2:length(lims)
        hs(1,k) = length(V(V>=lims(k-1)&V<lims(k))) / N;
    end
    hs(1,length(lims)+1) = length(V(V>=lims(end))) / N;
    speed_dist{sgs} = hs / sum(hs(:));
    
    edges = 1:0.1:5;
    nm = histc(Distance,edges);      % arrow length distribution
    nm = nm(1:end-1);
    nm_rows = size(Convergence,1);
    nm_cols = size(Convergence,2);
    chnum = nm_rows * nm_cols;
    dm = zeros(chnum,chnum);
    for x = 1:chnum
        for y = 1:chnum
            [i j] = gridind2sub(x,nm_cols);
            [k l] = gridind2sub(y,nm_cols);
            dm(x,y) = sqrt((i-k)^2 + (j-l)^2);
        end
    end
    nm2 = histc(dm(dm>0),edges);
    nm2 = nm2(1:end-1);
    nmm = nm ./ nm2;    % normalization by the abundance of a given distance
    nmm2 = nan2zero(nmm);
    normlength_dist{sgs} = nmm2 / sum(nmm2);
    length_dist{sgs} = nm / sum(nm);
    
    edges = 1.5:0.1:3.5;
    na = histc(Strength,edges);      % connection strength distribution
    na = na(1:end-1);
    strength_dist{sgs} = na / sum(na(:));
end

% Correlations - speed
next = 1;
D = zeros(1,3);
for xi = 1:2
    for yi = xi+1:3
        P = speed_dist{xi}(:);
        Q = speed_dist{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'speed.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,D','sheet1',['B' num2str(pref)])

% Correlations - arrow length
next = 1;
D = zeros(1,3);
for xi = 1:2
    for yi = xi+1:3
        P = length_dist{xi}(:);
        Q = length_dist{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'distance.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,D','sheet1',['B' num2str(pref)])

% Correlations - normalized arrow length
next = 1;
D = zeros(1,3);
for xi = 1:2
    for yi = xi+1:3
        P = normlength_dist{xi}(:);
        Q = normlength_dist{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'normdistance.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,D','sheet1',['B' num2str(pref)])

% Correlations - connection strength
next = 1;
D = zeros(1,3);
for xi = 1:2
    for yi = xi+1:3
        P = strength_dist{xi}(:);
        Q = strength_dist{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'strength.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,D','sheet1',['B' num2str(pref)])



% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end