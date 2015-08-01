function aphase_corr2
%APHASE_CORR2   Correlation between EEG frequency and firing pattern.
%   APHASE_CORR2 calculates and saves correlation between 
%       (i) mean phase angle of all first spikes and EEG frequency
%       (ii) burstiness and EEG frequency.
%
%   See also APHASE and APHASERUN_BURST2.

% Directories
global DATAPATH
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl_control\AUJ34_n2\';    % mat files
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\control\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\PhaseCorr2\'];
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running APHASE CORR2...');
global WB
WB(end+1) = wb;

% Main
sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
fre = [];
meanang = [];
brstness = [];
PPhist = [];
for o = 1:sf
    fname = files_short{o}     % filename
    cmps = strread(fname,'%s','delimiter','_');     % waitbar
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff = [inpdir1 fname];       % load
    load(ff)
    eeg = data(:,2)';
    len = length(data);
    clear data eeg0
    ff2 = [inpdir2 fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    seglen = 60 * sr;        % 60 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen -1;
    burstiness = zeros(1,lenr);
    burstfreq = zeros(1,lenr);
    frate = zeros(1,lenr);
    vburst = vdisc(Burst);
    Phist = zeros(length(edges)-1,lenr);
    meanangle = zeros(1,lenr);
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
        
        nm = histc(aang_afsp,edges);   % phase histogram
        Phist(:,k) = nm(1:end-1)';
        
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
    
    fre = [fre freq];
    meanang = [meanang meanangle];
    brstness = [brstness burstiness];
    PPhist = [PPhist Phist];
end
close(wb)

% Correlation coefficient
preR = corrcoef(meanang,fre);    % correlation between mean angle and EEG frequency
Rangfre = preR(2);
preR = corrcoef(brstness,fre);    % correlation between burstiness and EEG frequency
Rburstfre = preR(2);

% Plot
H1 = figure;
plot(meanang,fre,'.')
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.75+y_lim(1),num2str(Rangfre))
xlabel('mean angle')
ylabel('EEG frequency')
titlestr = [cmps{1} ' ' cmps{2}];
title(titlestr)
H2 = figure;
plot(brstness,fre,'.')
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.75+y_lim(1),num2str(Rburstfre))
xlabel('burstiness')
ylabel('EEG frequency')
title(titlestr)
H3 = figure;
subplot(4,1,1)
imagesc(flipud(PPhist))
title(titlestr)
subplot(4,1,2)
plot(meanang)
xlim([0,length(meanang)])
title('mean angle')
subplot(4,1,3)
plot(fre)
xlim([0,length(fre)])
title('EEG frequency')
subplot(4,1,4)
plot(brstness)
xlim([0,length(brstness)])
title('burstiness')

% Save
cd(resdir)
fns = [cmps{1} '_' cmps{2} '_ANGFRECORR.fig'];
saveas(H1,fns)
fns = [cmps{1} '_' cmps{2} '_BURSTFRECORR.fig'];
saveas(H2,fns)
fns = [cmps{1} '_' cmps{2} '_PHASECORR.fig'];
saveas(H3,fns)
fns = [cmps{1} '_' cmps{2} '_ANGFRECORR.mat'];
save(fns,'meanang','fre','Rangfre')
fns = [cmps{1} '_' cmps{2} '_BURSTFRECORR.mat'];
save(fns,'brstness','fre','Rburstfre')
close all

cd(mm)

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