function rjeta
%RJETA   EEG-trigerred spike train average.
%   RJETA calulates spike train overlays and firing histograms triggered by
%   negative peaks of the EEG.
%
%   See also RJRIPPLE2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\'];
inpdir2 = [DATADIR 'raphe_matfiles\raphe_juxta_files\'];
thetadir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\theta_segments\'];
nothdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\nontheta_segments\'];
sharpdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\raphe_juxta\STA_test\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RJETA...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-6) '.mat'];
    load(ff2)
    unit = Unit.values;
    eeg_times = linspace(1,length(eeg)/sr,length(eeg));
    
    bg = ldisc(eeg,sr);    % find negative peaks of the EEG
    ff = [nothdir 'NONTHETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % load non-theta segments
    load(ff)
    inxs = [];
    for k = 1:length(bg)    % restrict to non-theta segments
        [x y] = find(NonThetaSegments<bg(k),1,'last');
        if x == 2
            inxs = [inxs k];
        end
    end
    bg(inxs) = [];
    figure;plot(eeg)
    for k = 1:length(bg)
        line([bg(k) bg(k)],ylim,'Color','r')
    end
    line([0 length(eeg)],[-0.5 -0.5],'Color','g')
    
    if isempty(bg)
        continue
    end
    spwno = length(bg);
    pks = zeros(1,spwno);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for k = 1:spwno
        cnt = bg(k);
        if cnt - wn < 1 || cnt + wn > length(unit)
            continue
        end
        locunit = unit(cnt-wn:cnt+wn);
        S1 = subplot(2,1,1);
        hold on
        plot(locunit)
        locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
        z(locvdisc) = 1;
    end
    z2 = reshape(z,2*wn/20,20);
    subplot(2,1,2)
    bar(1:20,sum(z2))
    fn = [fname(1:end-4) '_ETA'];
    saveas(gcf,fn);
    
    close all
    waitbar(o/sf)
end
close(wb)
cd(mm)

% -------------------------------------------------------------------------
function vdisc = ldisc(eeg,sr)

% Threshold
figure
plot(eeg,'m');
title('Give the threshold! /Command window/')
thres = input('Give the threshold! ');

% Discriminating
unit = -eeg;
thres = -thres;
disc0 = find(unit>=thres);
segs = segm(disc0);
segs2 = uniteseg(segs,sr);
lvd = size(segs2,2);
vdisc = zeros(1,lvd);
for k = 1:lvd
    [maxe maxh] = max(unit(segs2(1,k):segs2(2,k)));
    vdisc(k) = segs2(1,k) + maxh - 1;
end


% discl = length(disc0);
% disc0 = [disc0; unit(disc0(1:discl))];
% disc0(1,discl+1) = length(unit) + 2;
% dif = diff(disc0(1,:));
% difn1 = find(dif>1);
% difn1(2:length(difn1)+1) = difn1;
% difn1(1) = 0;
% vdisc = zeros(1,length(difn1)-1);
% for j = 1:length(difn1) - 1
%     if difn1(j+1) - difn1(j) < 10
%         continue
%     end
%     [maxe,maxh] = max(disc0(2,difn1(j)+1:difn1(j+1)));
%     vdisc(j) = disc0(1,difn1(j)+maxh);
% end

% ----------------------------------------------------------------------------------------
function th = segm(ip)
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
preth = zeros(2,lenfdr+1);
preth(1,1) = ip(1);
for t = 1:lenfdr
    preth(2,t) = ip(fdr(t));
    preth(1,t+1) = ip(fdr(t)+1);
end
preth(2,end) = ip(end);

lpth = preth(2,:) - preth(1,:);
inx = find(lpth>10);
th = [preth(1,inx); preth(2,inx)];

% -------------------------------------------------------------------------
function segments2 = uniteseg(segments,sr)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    if (la-fi)/sr < 0.2    % eliminate <200 ms gaps
        [fnx fny] = find(segments2==fi);
        segments2(fnx,fny) = segments(2,k+1);
        segments2 = [segments2(1:2,1:fny) segments2(1:2,fny+2:end)];
    end
end