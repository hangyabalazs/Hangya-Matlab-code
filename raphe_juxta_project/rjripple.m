function rjripple
%RJRIPPLE   Firing of MRN cells during ripples.
%   RJRIPPLE calulates firing histograms, raw unit overlays and filtered
%   (90-140 Hz) LFP overlays centered on riplle maxima.
%
%   See also RAPHEVIEW2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\'];
inpdir2 = [DATADIR 'raphe_matfiles\raphe_juxta_files\'];
thetadir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\theta_segments\'];
nothdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\nontheta_segments\'];
sharpdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\raphe_juxta\Ripple\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RJRIPPLE...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
for o = 3:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-6) '.mat'];
    load(ff2)
    unit = Unit.values;
    
    ff = [sharpdir 'SHARPWAVE_SEGMENTS_' cmps{1} '_' ...
        cmps{2}(1:2)];  % sharpwaves; load sharpwave segments
    load(ff)
    flt = fir1(2048,[90/nqf 140/nqf]);
    feeg = filtfilt(flt,1,eeg);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for k = 1:size(SharpWaveSegments,2)
        big = SharpWaveSegments(1,k);
        fin = SharpWaveSegments(2,k);
        if fin > length(feeg)
            continue
        end
        lfeeg = feeg(big:fin);
        pk = find(lfeeg==max(lfeeg));
        cnt = big + pk(1);
        if cnt - wn < 1 || cnt + wn > length(feeg)
            continue
        end
        loceeg = feeg(cnt-wn:cnt+wn);
        subplot(3,1,1)
        hold on
        plot(loceeg)
        subplot(3,1,2)
        hold on
        locunit = unit(cnt-wn:cnt+wn);
        plot(locunit)
        locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
        z(locvdisc) = 1;
    end
    z2 = reshape(z,2*wn/20,20);
    subplot(3,1,3)
    bar(1:20,sum(z2))
    title(titlestr)
    fn = [fname(1:end-4) '_RIPPLE'];
    saveas(gcf,fn);
    
    close all
    waitbar(o/sf)
end
close(wb)
cd(mm)

% -------------------------------------------------------------------------
function segments2 = uniteseg(segments,sr)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    if (la-fi)/sr < 0.5
        [fnx fny] = find(segments2==fi);
        segments2(fnx,fny) = segments(2,k+1);
        segments2 = [segments2(1:2,1:fny) segments2(1:2,fny+2:end)];
    end
end



% ----------------------------------------------------------------------------------
function segments = short_killer(segments)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<30000);         % leaving segments shorter than 3 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];



% -------------------------------------------------------------------------
function sacr = racorr(ncc,bno)
%RACORR   Autocorrelation.
%   RACORR2(VD,BNO) calculates autocorrelogram for discriminated unit VD, 
%   using a +-1000 ms time window and BNO bins.

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;

% Autocorrelogram
zunit1 = zeros(1,length(round(nc))+5);
zunit1(round(nc)) = 1;
acr = xcorr(zunit1,1*sr);
acr(length(acr)/2+0.5) = [];
acr = reshape(acr,length(acr)/bno,bno);     % window: -200 ms - 200 ms
sacr = sum(acr);



% -------------------------------------------------------------------------
function H = stafig(sta,sta_index1,sta_index2,nn,wn,sr,titlestr)

time = linspace(-wn/sr/2,wn/sr/2,length(sta));
H = figure;
plot(time,sta,'LineWidth',1.5)
ach = allchild(H);     % figure title
ax = findobj(ach,'type','axes');
title(ax(end),titlestr)
x_lim = xlim;
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')