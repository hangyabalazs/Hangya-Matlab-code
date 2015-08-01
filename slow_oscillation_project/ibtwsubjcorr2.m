function ibtwsubjcorr2(pno1,pt1,egs1,pno2,pt2,egs2)
%IBTWSUBJCORR2   Between-subject correlations.
%   IBTWSUBJCORR2 calculates between-subject correlation of propagation
%   speed, distance and strength. See IMISORUN for details on the above
%   variables. Distribution distances are measured by symmetric
%   KL-divergence.
%
%   IBTWSUBJCORR2(PNO1,PT1,EGS1,PNO2,PT2,EGS2) requires 6 string input 
%   arguments:
%       PNO1, PNO2: patient IDs
%       PT1, PT2: patient names
%       EGS1, EGS2: EEG file IDs
%
%   See also IMISORUN and IMISORUN_MTX.

% Directories
global DATADIR
global DATAPATH
% pno = num2str(39);
% % pno = 'n1';
% pt = 'gaal'
% egs = [{'10'} {'14'} {'40'}];
pat1 = ['oiti' pno1 '_' pt1];
pat2 = ['oiti' pno2 '_' pt2];
resdir = [DATAPATH 'Ulbert\Summary\Btwcorr\'];

% Load
speed_dist1 = cell(1,3);
strength_dist1 = cell(1,3);
length_dist1 = cell(1,3);
normlength_dist1 = cell(1,3);
for sgs = 1:3
    eg = egs1{sgs};
    inpdir = [DATAPATH 'Ulbert\OITI_' pno1 '_EEG_' eg '\Mtx\'];
    
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
    speed_dist1{sgs} = hs / sum(hs(:));
    
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
    normlength_dist1{sgs} = nmm2 / sum(nmm2);
    length_dist1{sgs} = nm / sum(nm);
    
    edges = 1.5:0.1:3.5;
    na = histc(Strength,edges);      % connection strength distribution
    na = na(1:end-1);
    strength_dist1{sgs} = na / sum(na(:));
end

speed_dist2 = cell(1,3);
strength_dist2 = cell(1,3);
length_dist2 = cell(1,3);
normlength_dist2 = cell(1,3);
for sgs = 1:3
    eg = egs2{sgs};
    inpdir = [DATAPATH 'Ulbert\OITI_' pno2 '_EEG_' eg '\Mtx\'];
    
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
    speed_dist2{sgs} = hs / sum(hs(:));
    
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
    normlength_dist2{sgs} = nmm2 / sum(nmm2);
    length_dist2{sgs} = nm / sum(nm);
    
    edges = 1.5:0.1:3.5;
    na = histc(Strength,edges);      % connection strength distribution
    na = na(1:end-1);
    strength_dist2{sgs} = na / sum(na(:));
end

% Correlations - speed
next = 1;
D = zeros(1,3);
for xi = 1:3
    for yi = 1:3
        P = speed_dist1{xi}(:);
        Q = speed_dist2{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'bspeed.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat1},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,{pat2},'sheet1',['B' num2str(pref)])
xlswrite(xlsname,D','sheet1',['C' num2str(pref)])

% Correlations - arrow length
next = 1;
D = zeros(1,3);
for xi = 1:3
    for yi = 1:3
        P = length_dist1{xi}(:);
        Q = length_dist2{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'bdistance.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat1},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,{pat2},'sheet1',['B' num2str(pref)])
xlswrite(xlsname,D','sheet1',['C' num2str(pref)])

% Correlations - normalized arrow length
next = 1;
D = zeros(1,3);
for xi = 1:3
    for yi = 1:3
        P = normlength_dist1{xi}(:);
        Q = normlength_dist2{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'bnormdistance.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat1},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,{pat2},'sheet1',['B' num2str(pref)])
xlswrite(xlsname,D','sheet1',['C' num2str(pref)])

% Correlations - connection strength
next = 1;
D = zeros(1,3);
for xi = 1:3
    for yi = 1:3
        P = strength_dist1{xi}(:);
        Q = strength_dist2{yi}(:);
        D(next) = (KLdist(P,Q) + KLdist(Q,P)) / 2;
        next = next + 1;
    end
end
xlsname = [resdir 'bstrength.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{pat1},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,{pat2},'sheet1',['B' num2str(pref)])
xlswrite(xlsname,D','sheet1',['C' num2str(pref)])



% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end