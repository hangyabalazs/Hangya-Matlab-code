function noseg8
%NOSEG8   Analysis of theta and non-theta segments (NO project).
%   NOSEG8 calculates mean theta cycle length and its SD. The results are
%   plotted against time.
%
%   See also NOSEG7.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
inpdir2 = [DATAPATH 'NO\Wavelet\Segments\'];
resdir = [DATAPATH 'NO\Seg\'];
tblfile = [DATAPATH 'NO\seg_data.xls'];

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSEG8...','Position',[360 250 275 50]);
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
    
    ntf = ThetaSegments(1,:);   % theta cycle length and SD
    ntl = ThetaSegments(2,:);
    ntlen = (ntl - ntf) / sr;
    lenin = length(ntf);
    cyclen = zeros(1,lenin);    % theta cycle length
    cycsd = zeros(1,lenin);     % theta cycle length SD
    tms = zeros(1,lenin);
    for k = 1:lenin
        leeg = eeg(ntf(k):ntl(k));
        flt = fir1(256,[2.5 6]/nqf,'bandpass');      % bandpass filtering
        feeg = filtfilt(flt,1,leeg);
        feeg = (feeg - mean(feeg)) / std(feeg);
        ahee = angle(hilbert(feeg));    % Hilbert-transformation
        [cyclen(k) cycsd(k)] = eegfre(feeg,ahee,dsr);    % cycle length
        
        thinx1 = round(ntf(k)/newstep);
        thinx2 = round(ntl(k)/newstep);
        tms(k) = mean(thinx1,thinx2);
    end
    H1 = figure;    % plot
    title([fname(1:4) ' ' id2])
    plot(tms,cyclen,'.r','MarkerSize',16)
    H2 = figure;
    title([fname(1:4) ' ' id2])
    plot(tms,cycsd,'.r','MarkerSize',16)
    
    fntp = [resdir 'CycLen\' fname(1:end-4) '_MCYCLEN.fig'];    % save
    saveas(H1,fntp)
    close(H1)
    fntp = [resdir 'CycLen\' fname(1:end-4) '_SDCYCLEN.fig'];
    saveas(H2,fntp)
    close(H2)
    waitbar(o/sf)
end
close(wb)

% Save
% xlswrite([resdir 'noallpow_pre1800post3600.xls'],names','sheet1','A1')
% xlswrite([resdir 'noallpow_pre1800post3600.xls'],tps','sheet1','B1')

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
function [cyclen cycsd] = eegfre(feeg,ahee,sr)

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
cl6 = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    lg = (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.4 * sr);
    if ~lg
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
end
cyclen = mean(cl6);   % cycle length in ms
cycsd = std(cl6);   % SD of cycle length