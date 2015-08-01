function noseg3
%NOSEG3   Analysis of theta and non-theta segments (NO project).
%   NOSEG3 calculates the all spectral power (wavelet) of non-theta
%   segments in the last 1800 s for pre-injection files or in the first
%   3600 s in post-injection files. The results are plotted against time
%   and also saved in an Excel file.
%
%   See also NOSPECT.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
inpdir2 = [DATAPATH 'NO\Wavelet\Segments\'];
resdir = [DATAPATH 'NO\Seg\'];
tblfile = [DATAPATH 'NO\seg_data.xls'];
% mm = pwd;
% cd(resdir)

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSEG3...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Load segment data
[tbl0 tbl] = xlsread(tblfile);

% Main
sr = 5000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
newstep = 50;
dsr = sr / newstep;    % downsample on 100 Hz
names = cell(1,sf);
tps = zeros(1,sf);
for o = 1:sf
    fname = files(o).name;      % load data
    names{o} = fname;
    cmps = strread(fname(1:end-4),'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    
    load([inpdir fname])
    eeg = hEEG.values;
    len = length(eeg);
    
    fn2 = [fname(1:end-4) '_SEGMENTS'];     % load segment identification
    load([inpdir2 fn2])
    inx = find(strcmp({tbl{:,1}},fname));
    id1 = tbl{inx,2};
    id2 = tbl{inx,3};
    
    [poweeg,f] = eegwavelet(eeg(1:newstep:end),dsr);    % EEG wavelet
    pwind1_theta = find(f>6,1,'last');    % theta frequency band bounderies
    pwind2_theta = find(f<2.5,1,'first');
    pwind1_delta = find(f>2.5,1,'last');    % delta frequency band bounderies
    pwind2_delta = find(f<0.5,1,'first');
    thetapower = mean(poweeg(pwind1_theta:pwind2_theta,:));
    deltapower = mean(poweeg(pwind1_delta:pwind2_delta,:));
    allpower = mean(poweeg);
    
    ntf = NonThetaSegments(1,:);    % all power in non-theta segments
    ntl = NonThetaSegments(2,:);
    ntlen = (ntl - ntf) / sr;
    if strcmp(id2,'pre')
        ntinx = find(ntf>len-1800*sr);
    elseif strcmp(id2,'post')
        ntinx = find(ntl<3600*sr);
    end
    ntfi = ntf(ntinx);
    ntli = ntl(ntinx);
    lenin = length(ntinx);
    tp = zeros(1,lenin);
    H = figure;
    title([fname(1:4) ' ' id2])
    hold on
    for k = 1:lenin    % segment loop
        thinx1 = round(ntfi(k)/newstep);
        thinx2 = round(ntli(k)/newstep);
        tp(k) = mean(allpower(thinx1:thinx2));
        plot(mean(thinx1,thinx2),tp(k),'.r','MarkerSize',16)
    end
    tps(o) = mean(tp);
    fntp = [resdir 'AllPower\' fname(1:end-4) '_AP.fig'];
    saveas(H,fntp)      % save figure
    close(H)
    waitbar(o/sf)
end
close(wb)

% Save
xlswrite([resdir 'noallpow_pre1800post3600.xls'],names','sheet1','A1')
xlswrite([resdir 'noallpow_pre1800post3600.xls'],tps','sheet1','B1')

% cd(mm)

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [pow,f] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 1;
dj = 0.08;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0.4;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;