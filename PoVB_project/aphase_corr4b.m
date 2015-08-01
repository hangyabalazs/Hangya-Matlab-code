function aphase_corr4b(inpdir1)
%APHASE_CORR4B   Correlation between EEG frequency and firing pattern.
%   APHASE_CORR4B(DR) calculates and saves correlation between 
%       (i) mean phase angle of all first spikes and EEG frequency
%       (ii) burstiness and EEG frequency.
%   Input directory should be given as an argument (DR).
%
%   APHASE_CORR4B calculates 'bas' and 'bic' correlation for the same
%   cell.
%
%   See also APHASE and APHASERUN_BURST2.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
inpdir_bas = [inpdir1 'bas\'];
inpdir_bic = [inpdir1 'bic\'];
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\PhaseCorr4b_3\'];
mm = pwd;

% Filelist
[files1_bas files_short1_bas] = filelist(inpdir_bas);
[files1_bic files_short1_bic] = filelist(inpdir_bic);
[files2 files_short2] = filelist2(inpdir2);
files_short_bas = intersect(files_short1_bas,files_short2);
files_short_bic = intersect(files_short1_bic,files_short2);
sf_bas = length(files_short_bas);
sf_bic = length(files_short_bic);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running APHASE CORR4B...');
global WB
WB(end+1) = wb;

% Main
[meanang_bas, meanrl_bas, brstness_bas, fre_bas] = main(inpdir1,inpdir2,files_short_bas,sf_bas,wb,'bas');
[meanang_bic, meanrl_bic, brstness_bic, fre_bic, cmps] = main(inpdir1,inpdir2,files_short_bic,sf_bic,wb,'bic');

% Restrict frequency
fre_all = [fre_bas fre_bic];
sfr = sort(fre_all,'ascend');
lensfr = length(sfr);
pntno_bas = zeros(1,lensfr);
pntno_bic = zeros(1,lensfr);
pntno = zeros(1,lensfr);
c1 = 0.25;     % 0.25 Hz frequency band
for  k = 1:lensfr
    ind1 = sfr(k);
    ind2 = sfr(k) + c1;
    pntno_bas(k) = length(find(fre_bas>=ind1&fre_bas<=ind2));
    pntno_bic(k) = length(find(fre_bic>=ind1&fre_bic<=ind2));
    pntno(k) = min(pntno_bas(k),pntno_bic(k));
end
fn = find(pntno==max(pntno));
fn = fn(1);
ind1 = sfr(fn);
ind2 = sfr(fn) + c1;
lastp = sfr(find(sfr<ind2,1,'last'));
lag = ind2 - lastp;
fbind1 = ind1 - lag / 2;
fbind2 = ind2 - lag / 2;
pnts_bas = meanang_bas(fre_bas>fbind1&fre_bas<fbind2);
rad_bas = pnts_bas / 180 * pi;
pnts_bic = meanang_bic(fre_bic>fbind1&fre_bic<fbind2);
rad_bic = pnts_bic / 180 * pi;
fbftm_bas = sum(exp(1).^(i*rad_bas)) / length(rad_bas);    % first trigonometric moment
fbang_bas = angle(fbftm_bas) * 180 / pi;   % mean angle (deg.)
fbftm_bic = sum(exp(1).^(i*rad_bic)) / length(rad_bic);    % first trigonometric moment
fbang_bic = angle(fbftm_bic) * 180 / pi;   % mean angle (deg.)
try
    [U2 fbp] = b_watsontwo(rad_bas',rad_bic');
catch
    fbp = [NaN NaN];
end

% Restrict phase
ang_all = [meanang_bas meanang_bic];
sma = sort(ang_all,'ascend');
lensma = length(sma);
pntno_bas = zeros(1,lensma);
pntno_bic = zeros(1,lensma);
pntno = zeros(1,lensma);
c2 = 30;     % 30 deg. phase band
for  k = 1:lensma
    ind1 = sma(k);
    ind2 = sma(k) + c2;
    pntno_bas(k) = length(find(meanang_bas>=ind1&meanang_bas<=ind2));
    pntno_bic(k) = length(find(meanang_bic>=ind1&meanang_bic<=ind2));
    pntno(k) = min(pntno_bas(k),pntno_bic(k));
end
fn = find(pntno==max(pntno));
fn = fn(1);
ind1 = sma(fn);
ind2 = sma(fn) + c2;
lastp = sma(find(sma<ind2,1,'last'));
lag = ind2 - lastp;
pbind1 = ind1 - lag / 2;
pbind2 = ind2 - lag / 2;
pnts_bas = fre_bas(meanang_bas>pbind1&meanang_bas<pbind2);
pnts_bic = fre_bic(meanang_bic>pbind1&meanang_bic<pbind2);
pbfre_bas = mean(pnts_bas);   % mean frequency
pbfre_bic = mean(pnts_bic);
try
    [pbp,W] = ranksum(pnts_bas,pnts_bic);
catch
    pbp = NaN;
end

% Plot
close(wb)
H1 = figure;
plot(meanang_bas,fre_bas,'.')   % mean angle vs. EEG frequency
[gr icp err] = linefit(meanang_bas,fre_bas);
x = [min(meanang_bas):0.01:max(meanang_bas)];
y = x .* gr + icp;
hold on
plot(x,y)
plot(meanang_bic,fre_bic,'r.')
[gr icp err] = linefit(meanang_bic,fre_bic);
x = [min(meanang_bic):0.01:max(meanang_bic)];
y = x .* gr + icp;
plot(x,y,'r')
x_lim = xlim;
y_lim = ylim;
line([x_lim(1) x_lim(2)],[fbind1 fbind1],'Color','green')
line([x_lim(1) x_lim(2)],[fbind2 fbind2],'Color','green')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['bas mean: ' num2str(fbang_bas)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    ['bic mean: ' num2str(fbang_bic)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    ['p < ' num2str(fbp(2))],'Color','black')
line([pbind1 pbind1],[y_lim(1) y_lim(2)],'Color','green')
line([pbind2 pbind2],[y_lim(1) y_lim(2)],'Color','green')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.7,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['bas mean: ' num2str(pbfre_bas)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.7,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    ['bic mean: ' num2str(pbfre_bic)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.7,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    ['p = ' num2str(pbp)],'Color','black')
xlabel('mean angle')
ylabel('EEG frequency')
titlestr = [cmps{1} ' ' cmps{2}];
title(titlestr)

H2 = figure;
plot(brstness_bas,fre_bas,'.')      % burstiness vs. EEG frequency
% [gr icp err] = linefit(brstness_bas,fre_bas);
% x = [min(brstness_bas):0.01:max(brstness_bas)];
% y = x .* gr + icp;
hold on
% plot(x,y)
plot(brstness_bic,fre_bic,'r.')
% [gr icp err] = linefit(brstness_bic,fre_bic);
% x = [min(brstness_bic):0.01:max(brstness_bic)];
% y = x .* gr + icp;
% plot(x,y,'r')
axis([0 1 0.5 3])
xlabel('burstiness')
ylabel('EEG frequency')
title(titlestr)

H3 = figure;
plot(meanrl_bas,fre_bas,'.')      % mean res. length vs. EEG frequency
% [gr icp err] = linefit(meanrl_bas,fre_bas);
% x = [min(meanrl_bas):0.01:max(meanrl_bas)];
% y = x .* gr + icp;
hold on
% plot(x,y)
plot(meanrl_bic,fre_bic,'r.')
% [gr icp err] = linefit(meanrl_bic,fre_bic);
% x = [min(meanrl_bic):0.01:max(meanrl_bic)];
% y = x .* gr + icp;
% plot(x,y,'r')
axis([0 1 0.5 3])
xlabel('mean res. length')
ylabel('EEG frequency')
title(titlestr)

H4 = figure;
plot(meanang_bas,meanrl_bas,'.')      % mean phase vs. mean res. length
% [gr icp err] = linefit(meanang_bas,meanrl_bas);
% x = [min(meanang_bas):0.01:max(meanang_bas)];
% y = x .* gr + icp;
hold on
% plot(x,y)
plot(meanang_bic,meanrl_bic,'r.')
% [gr icp err] = linefit(meanang_bic,meanrl_bic);
% x = [min(meanang_bic):0.01:max(meanang_bic)];
% y = x .* gr + icp;
% plot(x,y,'r')
ylim([0 1])
xlabel('mean angle')
ylabel('mean res. length')
title(titlestr)

% Save
cd(resdir)
if length(cmps) > 2
    str = [cmps{1} '_' cmps{2} '_' cmps{3}];
else
    str = [cmps{1} '_' cmps{2}];
end
fns = [str '_ANGFRECORR.fig'];
saveas(H1,fns)
fns = [str '_ANGFRECORR.jpg'];
saveas(H1,fns)
fns = [str '_BURSTFRECORR.fig'];
saveas(H2,fns)
fns = [str '_BURSTFRECORR.jpg'];
saveas(H2,fns)
fns = [str '_MVLFRECORR.fig'];
saveas(H3,fns)
fns = [str '_MVLFRECORR.jpg'];
saveas(H3,fns)
fns = [str '_ANGMVLCORR.fig'];
saveas(H4,fns)
fns = [str '_ANGMVLCORR.jpg'];
saveas(H4,fns)
fns = [str '_ANGFRECORR.mat'];
save(fns,'meanang_bas','fre_bas','meanang_bic','fre_bic')
fns = [str '_BURSTFRECORR.mat'];
save(fns,'brstness_bas','fre_bas','brstness_bic','fre_bic')
fns = [str '_MVLFRECORR.mat'];
save(fns,'meanrl_bas','fre_bas','meanrl_bic','fre_bic')
fns = [str '_ANGMVLCORR.mat'];
save(fns,'meanang_bas','meanrl_bas','meanang_bic','meanrl_bic')
close all

cd(mm)

% -------------------------------------------------------------------------
function [meanang, meanrl, brstness, fre, cmps] = main(inpdir1,inpdir2,files_short,sf,wb,bob);

sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
fre = [];
meanang = [];
meanrl = [];
brstness = [];
PPhist = [];
for o = 1:sf
    fname = files_short{o}     % filename
    cmps = strread(fname(1:end-4),'%s','delimiter','_');     % waitbar
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff = [inpdir1 bob '\' fname];       % load
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
    frate = zeros(1,lenr);
    vburst = vdisc(Burst);
    Phist = zeros(length(edges)-1,lenr);
    meanangle = zeros(1,lenr);
    mrl = zeros(1,lenr);
    freq = zeros(1,lenr);
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);
        loceeg = eeg(ind1(k):ind2(k));
        
% Burst statistics
        lvb = locvburst(vburst,ind1(k),ind2(k));
        burstnum = size(lvb,2);
        intraburstiv = [];
        intraburstnum = zeros(1,burstnum);
        for j = 1:burstnum      % computing intraburstiv
            b = vdisc(vdisc>=lvb(1,j)&vdisc<=lvb(2,j));
            db = diff(b);
            intraburstiv = [intraburstiv db];
            intraburstnum(j) = length(b);   %intraburst spike number
        end
        burstiness(k) = (length(intraburstiv) + burstnum) / length(vd);
        burstlength = (lvb(2,:) - lvb(1,:)) / sr;
        if ~isempty(intraburstnum)
            intraburstfreq(k) = mean((intraburstnum-1)./burstlength);
            ibspno(k) = mean(intraburstnum);
        else
            intraburstfreq(k) = NaN;
            ibspno(k) = NaN;
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
        
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
    
    fre = [fre freq];
    meanang = [meanang meanangle];
    meanrl = [meanrl mrl];
    brstness = [brstness burstiness];
    PPhist = [PPhist Phist];
end

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist2(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-11) '.mat'];
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function vburst = locvburst(vburst,ind1,ind2)

VBurst = vburst;
vburst = VBurst .* (VBurst>ind1) .* (VBurst<ind2);
if isequal(size(VBurst),[1,2])
    VBurst = VBurst';
end
vburst = VBurst .* (VBurst>ind1) .* (VBurst<ind2);  % restrict vburst to the window
if isempty(vburst)
    return
end
vb1 = [vburst(1,:)>0 vburst(2,:)>0] .* [[1:size(vburst,2)] (-1)*[1:size(vburst,2)]];
if sum(vb1) > 0     % handle the case when burst exceeds the window
    ind = find(vburst,1,'last')+1;
    vburst(ind) = ind2;
    vb1 = [vburst(1,:)>0 vburst(2,:)>0] .* [[1:size(vburst,2)] (-1)*[1:size(vburst,2)]];
end
if sum(vb1) < 0
    ind = find(vburst,1,'first')-1;
    vburst(ind) = ind1;
end
fV1 = find(VBurst>ind2,1,'first');
fV2 = find(VBurst<ind1,1,'last');
if fV1 - fV2 == 1 & mod(fV1,2) == 0
    vburst = [ind1 ind2];   % when one burst exceeds the window on both sides
end
vburst = vburst(find(vburst));
vburst = reshape(vburst,2,length(vburst)/2);

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