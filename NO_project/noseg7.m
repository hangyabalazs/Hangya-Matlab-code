function noseg7
%NOSEG7   Analysis of theta and non-theta segments (NO project).
%   NOSEG7 calculates the frequency, amplitude (SD), mean and total
%   duration of non-theta, theta and transition segments in 15 min
%   overlapping windows. The results are plotted against time.
%
%   See also NOSEG and NOSEG2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
inpdir2 = [DATAPATH 'NO\Wavelet\Segments\'];
resdir1 = [DATAPATH 'NO\Seg2\NonThetaFreq\'];
resdir2 = [DATAPATH 'NO\Seg2\ThetaFreq\'];
resdir3 = [DATAPATH 'NO\Seg2\TransitionFreq\'];
resdir4 = [DATAPATH 'NO\Seg2\NonThetaPlusTransitionFreq\'];
tblfile = [DATAPATH 'NO\seg_data.xls'];

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSEG7...','Position',[360 250 275 50]);
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
for o = 2:sf
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
    
    ntf = NonThetaSegments(1,:);    % segment limits
    ntl = NonThetaSegments(2,:);
    ntlen = (ntl - ntf) / sr;
    tf = ThetaSegments(1,:);
    tl = ThetaSegments(2,:);
    tlen = (tl - tf) / sr;
    tsf = Transitions(1,:);
    tsl = Transitions(2,:);
    tslen = (tsl - tsf) / sr;
    
    winlen = 15 * 60 * sr;   % window size (15 min)
    [k1 k2] = size(eeg');
    maxi = floor(k2/winlen);
    ovlp = 15;
    mx = floor((k2-winlen)/(winlen/ovlp));
    ntsds = zeros(1,mx);
    tsds = zeros(1,mx);
    ntlens = zeros(1,mx);
    tlens = zeros(1,mx);
    tslens = zeros(1,mx);
    ntlenscum = zeros(1,mx);
    tlenscum = zeros(1,mx);
    tslenscum = zeros(1,mx);
    ntfreq = zeros(1,mx);
    ntfreq2 = zeros(1,mx);
    tfreq = zeros(1,mx);
    tsfreq = zeros(1,mx);
    tm = zeros(1,mx);
    for t = 1:mx
        inx1 = (t - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
        inx1 = round(inx1);
        inx2 = inx1 + winlen - 1;
        tm(t) = (inx1 + inx2) / 2;
        
        ntinx = find(ntf>inx1&ntl<inx2);
        ntfi = ntf(ntinx);
        ntli = ntl(ntinx);
        ntleni = ntlen(ntinx);
        lenin1 = length(ntinx);
        ntsd = zeros(1,lenin1);
        for k = 1:lenin1
            ntsd(k) = std(eeg(ntfi(k):ntli(k)));
        end
        ntsds(t) = mean(ntsd);  % mean non-theta amp. (SD)
        ntlens(t) = mean(ntleni);  % mean non-theta segment length
        ntlenscum(t) = sum(ntleni) / (winlen / sr);  % proportion of non-theta segments
        ntfreq(t) = length(ntinx) / (winlen / sr);  % frequency of non-theta segments
        
        tinx = find(tf>inx1&tl<inx2);
        tfi = tf(tinx);
        tli = tl(tinx);
        tleni = tlen(tinx);
        lenin2 = length(tinx);
        tsd = zeros(1,lenin2);
        for k = 1:lenin2
            tsd(k) = std(eeg(tfi(k):tli(k)));
        end
        tsds(t) = mean(tsd);  % mean theta amp. (SD)
        tlens(t) = mean(tleni);  % mean theta segment length
        tlenscum(t) = sum(tleni) / (winlen / sr);  % proportion of theta segments
        tfreq(t) = length(tinx) / (winlen / sr);  % frequency of theta segments
        
        tsinx = find(tsf>inx1&tsl<inx2);
        tsfi = tsf(tsinx);
        tsli = tsl(tsinx);
        tsleni = tslen(tsinx);
        tslens(t) = mean(tsleni);  % mean transition amp. (SD)
        tslenscum(t) = sum(tsleni) / (winlen / sr);  % mean transition segment length
        tsfreq(t) = length(tsinx) / (winlen / sr);  % proportion of transition segments
        ntfreq2(t) = (length(ntinx) + length(tsinx)) / (winlen / sr);  % frequency of transition segments
    end

% Plot and save
%     H = figure;
%     plot(tm/sr,ntsds./tsds)
%     ylim([0 6])
%     title(titlestr)
%     fnp = [resdir 'NonThetaPerThetaSD_' fname(1:end-4) '.fig'];
%     saveas(H,fnp)
    H = figure;
    plot(tm/sr,ntfreq)
    ylim([0 0.85])
    title(titlestr)
    fnp = [resdir1 'NonThetaFreq_' fname(1:end-4) '.fig'];
    saveas(H,fnp)
    H = figure;
    plot(tm/sr,ntfreq2)
    ylim([0 0.85])
    title(titlestr)
    fnp = [resdir4 'NonThetaPlusTransitionFreq_' fname(1:end-4) '.fig'];
    saveas(H,fnp)
    H = figure;
    plot(tm/sr,tfreq)
    ylim([0 0.3])
    title(titlestr)
    fnp = [resdir2 'ThetaFreq_' fname(1:end-4) '.fig'];
    saveas(H,fnp)
    H = figure;
    plot(tm/sr,tsfreq)
    ylim([0 0.6])
    title(titlestr)
    fnp = [resdir3 'TransitionFreq_' fname(1:end-4) '.fig'];
    saveas(H,fnp)
end
close(wb)

% Save
% xlswrite([resdir 'nonthetasd.xls'],names','sheet1','A1')
% xlswrite([resdir 'nonthetasd.xls'],ntsds','sheet1','B1')

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