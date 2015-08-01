function afre_corr(inpdir1)
%AFRE_CORR   Correlation between EEG frequency and firing pattern.
%   AFRE_CORR(DR) calculates and saves correlation between EEG frequency
%       (i) mean phase angle of all first spikes (circular-linear corr.)
%       (ii) mean resultant length
%       (iii) STA index
%       (iv) burstiness
%       (v) intraburst spike number
%       (vi) intraburst frequency.
%   Scatter-plots with regression F-test and Monte Carlo permutation test 
%   results, time-domain graphs and mat files containing the correlated
%   variables are saved.
%
%   Input directory should be given as an argument (DR). Edit program code
%   to modify result directory!
%
%   See also APHASERUN_BURST2_STAND and AFRECORR_CALL.

% Input argument check
error(nargchk(1,1,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\control\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\AfreCorr3\'];
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running AFRE CORR...');
global WB
WB(end+1) = wb;

% Main
[meanang, meanrl, PPhist, STAindex2, SSta, brstness, ibspikeno, ibfr,...
    firstisirec, fre cmps] = main(inpdir1,inpdir2,files_short,sf,wb);

% Plot
close(wb)
titlestr = [cmps{1} ' ' cmps{2}];
meanang_rad = meanang / 180 * pi;
H1 = circscplot(fre,meanang_rad,'EEG frequency','mean angle',titlestr);
H2 = linscplot(fre,meanrl,'EEG frequency','mean resultant length',titlestr);
H3 = linscplot(fre,STAindex2,'EEG frequency','STA index',titlestr);
H4 = linscplot(fre,brstness,'EEG frequency','burstiness',titlestr);
H5 = linscplot(fre,ibspikeno,'EEG frequency','intraburst spike number',titlestr);
H6 = linscplot(fre,ibfr,'EEG frequency','intraburst frequency',titlestr);
H61 = linscplot(fre,firstisirec,'EEG frequency',...
    'reciprocal of first intraburst interval',titlestr);

H7 = figure;
subplot(4,1,1)
plot(fre)
xlim([0,length(fre)])
title([titlestr ' EEG frequency'])
subplot(4,1,2)
imagesc(flipud(PPhist))
title('phase histogram')
subplot(4,1,3)
plot(meanang)
xlim([0,length(meanang)])
title('mean angle')
subplot(4,1,4)
plot(meanrl)
xlim([0,length(meanrl)])
title('mean resultant length')

H8 = figure;
subplot(3,1,1)
plot(fre)
xlim([0,length(fre)])
title([titlestr ' EEG frequency'])
subplot(3,1,2)
imagesc(flipud(SSta))
title('STA')
subplot(3,1,3)
plot(STAindex2)
xlim([0,length(STAindex2)])
title('STA index')

H9 = figure;
subplot(4,1,1)
plot(fre)
xlim([0,length(fre)])
title([titlestr ' EEG frequency'])
subplot(4,1,2)
plot(brstness)
xlim([0,length(brstness)])
title('burstiness')
subplot(4,1,3)
plot(ibspikeno)
xlim([0,length(ibspikeno)])
title('intraburst spike number')
subplot(4,1,4)
plot(ibfr)
xlim([0,length(ibfr)])
title('intraburst frequency')

% Save
cd(resdir)
dbclear if error
if length(cmps) > 2
    str = [cmps{1} '_' cmps{2} '_' cmps{3}];
else
    str = [cmps{1} '_' cmps{2}];
end
fns = [str '_FREANGCORR.fig'];
saveas(H1,fns)
fns = [str '_FREANGCORR.eps'];
saveas(H1,fns)
fns = [str '_FREMVLCORR.fig'];
saveas(H2,fns)
fns = [str '_FREMVLCORR.eps'];
saveas(H2,fns)
fns = [str '_FRESTACORR.fig'];
saveas(H3,fns)
fns = [str '_FRESTACORR.eps'];
saveas(H3,fns)
fns = [str '_FREBURSTINESSCORR.fig'];
saveas(H4,fns)
fns = [str '_FREBURSTINESSCORR.eps'];
saveas(H4,fns)
fns = [str '_FREIBSPNOCORR.fig'];
saveas(H5,fns)
fns = [str '_FREIBSPNOCORR.eps'];
saveas(H5,fns)
fns = [str '_FREIBFRCORR.fig'];
saveas(H6,fns)
fns = [str '_FREIBFRCORR.eps'];
saveas(H6,fns)
fns = [str '_FREFISIRCORR.fig'];
saveas(H61,fns)
fns = [str '_FREFISIRCORR.eps'];
saveas(H61,fns)
fns = [str '_FREANGMVL.fig'];
saveas(H7,fns)
fns = [str '_FREANGMVL.eps'];
saveas(H7,fns)
fns = [str '_FRESTASTAINX.fig'];
saveas(H8,fns)
fns = [str '_FRESTASTAINX.eps'];
saveas(H8,fns)
fns = [str '_FREBURSTPARAM.fig'];
saveas(H9,fns)
fns = [str '_FREBURSTPARAM.eps'];
saveas(H9,fns)

fns = [str '_FREANGCORR.mat'];
save(fns,'meanang','meanrl','PPhist','STAindex2','SSta','fre')
fns = [str '_FREBURSTCORR.mat'];
save(fns,'brstness','ibspikeno','ibfr','fre')

close all
cd(mm)

% -------------------------------------------------------------------------
function [meanang, meanrl, PPhist, STAindex2, SSta, brstness, ibspikeno, ...
    ibfr, firstisirec, fre, cmps] = main(inpdir1,inpdir2,files_short,sf,wb)

sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
fre = [];
meanang = [];
meanrl = [];
brstness = [];
PPhist = [];
ibspikeno = [];
ibfr = [];
firstisirec = [];
STAindex2 = [];
SSta = [];
wn = 2 * dsr;    % 2 sec. window for STA
for o = 1:sf
    fname = files_short{o}     % filename
    cmps = strread(fname(1:end-4),'%s','delimiter','_');     % waitbar
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff = [inpdir1 '\' fname];       % load
    load(ff)
    eeg = data(:,2)';
    len = length(data);
    clear data eeg0
    ff2 = [inpdir2 fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    seglen = 30 * sr;        % 30 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen -1;
    burstiness = zeros(1,lenr);
    burstfreq = zeros(1,lenr);
    intraburstfreq = zeros(1,lenr);
    ibspno = zeros(1,lenr);
    firstisireciprocal = zeros(1,lenr);
    frate = zeros(1,lenr);
    vburst = vdisc(Burst);
    Phist = zeros(length(edges)-1,lenr);
    meanangle = zeros(1,lenr);
    mrl = zeros(1,lenr);
    sta_index1 = zeros(1,lenr);
    sta_index2 = zeros(1,lenr);
    Sta = zeros(wn+1,lenr);
    freq = zeros(1,lenr);
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);
        loceeg = eeg(ind1(k):ind2(k));
        
% Burst statistics
        lvb = locvburst(vburst,ind1(k),ind2(k));
        burstnum = size(lvb,2);
        intraburstiv = [];
        intraburstnum = zeros(1,burstnum);
        firstintraburstiv = zeros(1,burstnum);
        for j = 1:burstnum      % computing intraburstiv
            b = vdisc(vdisc>=lvb(1,j)&vdisc<=lvb(2,j));
            db = diff(b);
            intraburstiv = [intraburstiv db];
            firstintraburstiv(j) = db(1);
            intraburstnum(j) = length(b);   %intraburst spike number
        end
        burstiness(k) = (length(intraburstiv) + burstnum) / length(vd);
        burstlength = (lvb(2,:) - lvb(1,:)) / sr;
        if ~isempty(intraburstnum)
            intraburstfreq(k) = mean((intraburstnum-1)./burstlength);
            ibspno(k) = mean(intraburstnum);
            firstisireciprocal(k) = mean(1./firstintraburstiv);
        else
            intraburstfreq(k) = NaN;
            ibspno(k) = NaN;
            firstisireciprocal(k) = NaN;
        end
        burstfreq(k) = sr * (burstnum - 1)  / (vdisc(Burst(2,end)) - vdisc(Burst(1,1)));
        efflen = (vd(end) - vd(1)) / sr;
        frate(k) = (length(vd) - 1) / efflen;
        
% Phase histograms
        eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/const);
        
        ssi = vdisc;       % allfirstspikes
        for k2 = 1:size(Burst,2)
            ssi(Burst(1,k2)+1:Burst(2,k2)) = 0;
        end
        afsp = ssi(ssi>0);
        lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
        lafsp = round(lafsp/const);    % downsample unit on 1000 Hz
    
        [aang_afsp dinx_afsp cyclen_afsp] = laphase(eeg2,lafsp,dsr);    % PHASE - all first spikes
        n_afsp = length(aang_afsp);
        ftm_afsp = sum(exp(1).^(i*aang_afsp)) / n_afsp;    % first trigonometric moment
        ang_afsp = angle(ftm_afsp);   % mean angle
        mvl_afsp = abs(ftm_afsp);     % mean resultant length
        aang_afsp = aang_afsp * 180 / pi;
        ang_afsp = ang_afsp * 180 / pi;
        freq(k) = 1 / cyclen_afsp * 1000;
        meanangle(k) = ang_afsp;
        mrl(k) = mvl_afsp;
        
        nm = histc(aang_afsp,edges);   % phase histogram
        Phist(:,k) = nm(1:end-1)';
        
% Spike Triggered Average
        [sta sta_index1(k) sta_index2(k)] = astanorm(vdisc2,eeg2,wn);    % STA
        Sta(:,k) = sta';
        
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
    
    fre = [fre freq];
    meanang = [meanang meanangle];
    meanrl = [meanrl mrl];
    brstness = [brstness burstiness];
    ibspikeno = [ibspikeno ibspno];
    ibfr = [ibfr intraburstfreq];
    firstisirec = [firstisirec firstisireciprocal];
    PPhist = [PPhist Phist];
    STAindex2 = [STAindex2 sta_index2];
    SSta = [SSta Sta];
end

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct([]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        if isempty(files2)
            files2 = files(i);
        else
            files2(end+1) = files(i);
        end
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist2(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct([]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        if isempty(files2)
            files2 = files(i);
        else
            files2(end+1) = files(i);
        end
        files2_short{end+1} = [files(i).name(1:end-11) '.mat'];
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function lvburst = locvburst(vburst,ind1,ind2)

fst = vburst(1,:);
lst = vburst(2,:);
fi = find(fst>ind1,1,'first');
li = find(lst<ind2,1,'last');
lvburst = vburst(:,fi:li);

% -------------------------------------------------------------------------
function [ang inx cyclen] = laphase(eeg,vdisc,sr)
%APHASE    Phase angles for unit relative to EEG.
%   [A I C] = APHASE2(EEG,VDISC,SR) calculates Hilbert phase angles (A) for 
%   discriminated unit (VDISC) relative to EEG, when sampling frequency
%   is given in SR. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 100 ms length (half
%   wavelength of filter cutoff freq.). Indices of disclosed spikes of
%   vdisc are returned in I and mean cycle length in C.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 100 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
inx = find(vdisc<fn(1));
cl6 = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    else
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);
cyclen = mean(cl6) / sr * 1000;   % cycle length in ms

% -------------------------------------------------------------------------
function H = linscplot(x,y,xlb,ylb,titlestr)

H = figure;
plot(x,y,'.')
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
p2 = linrandomisationtest(x,y,'twosided');   % Monte Carlo randomisation test for sign. corr.

x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.6+y_lim(1),['p (perm. test): ' num2str(p2)])
xlabel(xlb)
ylabel(ylb)
title(titlestr)

% -------------------------------------------------------------------------
function H = circscplot(x,y,xlb,ylb,titlestr)

H = figure;
plot(x,y/pi*180,'.')
hold on
plot(x,y/pi*180+360,'.')
[b,bint,r,rint,stats] = regress(x',[ones(length(y),1) sin(y') cos(y')]);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
R = lincirc_corr(x,y);
R - R1
p2 = circrandomisationtest(x,y,'twosided');   % Monte Carlo randomisation test for sign. corr.

x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.6+y_lim(1),['p (perm. test): ' num2str(p2)])
xlabel(xlb)
ylabel(ylb)
title(titlestr)

% -------------------------------------------------------------------------
function T = linteststat(x,y,mnx,mny,n,str)
% Test statistic.

if isequal(str,'onesided')
    T = sum(x.*y) - n * mnx * mny;
else
    T = abs(sum(x.*y)-n*mnx*mny);
end

% -------------------------------------------------------------------------
function p = linrandomisationtest(x,y,str)

n1 = length(x);
mnx = mean(x);
mny = mean(y);
mno = 10000;     % maximal number of test statistics calculation
nk = factorial(n1);
if nk < mno
    N = nk;
    allcomb = y(perms(1:n1));     % all possible permutations
else
    N = mno;
    allcomb = zeros(mno,n1);
    rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
    for k = 1:mno
        allcomb(k,:) = y(randperm(n1));
    end
    disp('Randomisation distribution is estimated.')
end
frs = zeros(1,N);
for k = 1:N
    tng1 = x;
    tng2 = allcomb(k,:);
    frs(k) = linteststat(tng1,tng2,mnx,mny,n1,str);
end
fr = linteststat(x,y,mnx,mny,n1,str);
sfrs = sort(frs,'ascend');
if fr > sfrs(end)
    m = mno;
else
    m2 = find(sfrs>fr,1,'first');
    m1 = m2 - 1;
    m = m1 + (fr - sfrs(m1)) / (sfrs(m2) - sfrs(m1));
end
p = (N - m + 1) / N;

% -------------------------------------------------------------------------
function T = circteststat(x,y,str)
% Test statistic.

if isequal(str,'onesided')
    T = lincirc_corr(x,y);
else
    T = abs(lincirc_corr(x,y));
end

% -------------------------------------------------------------------------
function p = circrandomisationtest(x,y,str)

n1 = length(x);
mno = 10000;     % maximal number of test statistics calculation
nk = factorial(n1);
if nk < mno
    N = nk;
    allcomb = y(perms(1:n1));     % all possible permutations
else
    N = mno;
    allcomb = zeros(mno,n1);
    rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
    for k = 1:mno
        allcomb(k,:) = y(randperm(n1));
    end
    disp('Randomisation distribution is estimated.')
end
frs = zeros(1,N);
for k = 1:N
    tng1 = x;
    tng2 = allcomb(k,:);
    frs(k) = circteststat(tng1,tng2,str);
end
fr = circteststat(x,y,str);
sfrs = sort(frs,'ascend');
if fr > sfrs(end)
    m = mno;
else
    m2 = find(sfrs>fr,1,'first');
    m1 = m2 - 1;
    m = m1 + (fr - sfrs(m1)) / (sfrs(m2) - sfrs(m1));
end
p = (N - m + 1) / N;