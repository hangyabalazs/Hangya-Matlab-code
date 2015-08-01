function L5_spiking_5percent(inpdir1)
%AFRE_RESTRICT_STAND_LAYER5   Phase and burst analysis restricted to a frequency band.
%   AFRE_RESTRICT_STAND_LAYER5(DR) calculates and saves phase and burst
%   analysis results for a given EEG frequency band (1.1 - 1.6 Hz). See
%   ACLUSTERCUT and APHASE_BURST2_STAND for details on the analysis. Input
%   directory should be given as an argument (DR). EEG is standardized for
%   phase calculations.
%
%   See also ACLUSTERCUT and APHASERUN_BURST2_STAND.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
resdir1 = [DATAPATH 'Hajni_layer_5\L5_spiking_5percent\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir1);
sf = length(files_short);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running L5 SPIKING 5PERCENT...');
global WB
WB(end+1) = wb;

% Main
main(inpdir1,resdir1,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir1,resdir1,files_short,sf,wb)

sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
aang_fs = [];
aang_sp = [];
aang_as = [];
aang_afsp = [];
aang_ibsang = [];
aang_sspoang = [];
aang_allang = [];
cycnb = [];
r_fs = [];
r_sp = [];
r_as = [];
r_afsp = [];
H1 = figure;
for o = 1:sf
    fname = files_short{o}     % filename
    cmps = strread(fname(1:end-4),'%s','delimiter','_');     % waitbar
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff = [inpdir1 fname];       % load
    load(ff)
    len = length(eeg);
    
    nqf = dsr / 2;      % filtering EEG
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,eeg(1:const:end));
    feeg = (feeg - mean(feeg)) / std(feeg);
    ahee = angle(hilbert(feeg));    % Hilbert-transformation
    
    fst = Burst(1,:);   % burst first spikes
    vb1 = vdisc(fst);
    sso = vdisc;       % single spikes
    ssi = vdisc;       % allfirstspikes
    ib = zeros(size(vdisc));       % intraburst spikes
    for k = 1:size(Burst,2)
        sso(Burst(1,k):Burst(2,k)) = 0;
        ssi(Burst(1,k)+1:Burst(2,k)) = 0;
        ib(Burst(1,k):Burst(2,k)) = 1;
    end
    sspo = sso(sso>0);
    afsp = ssi(ssi>0);
    ibs = vdisc(find(ib));
    vburst = vdisc(Burst);
    
    seglen = 30 * sr;        % 30 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen -1;
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);      % localize
        loceeg = eeg(ind1(k):ind2(k));
        lfeeg = feeg((ind1(k)-1)/const+1:ind2(k)/const);
        lahee = ahee((ind1(k)-1)/const+1:ind2(k)/const);
        
% Phase histograms
        eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/const);
        
        cyclen = eegfre(lfeeg,lahee,dsr);    % filter EEG, Hilbert-transform
        freq = 1 / cyclen * 1000;
        frlim1 = 1.1
        frlim2 = 1.6
        if frlim1 < freq & freq < frlim2        % restrict frequency band
            lvb1 = vb1(vb1>ind1(k)&vb1<ind2(k)) - ind1(k);
            lvb1 = round(lvb1/const);    % burst first spikes, downsample unit on 1000 Hz
            figure
            plot(eeg2,'g')
            hold on
            [paang_fs pdinx_fs] = laphase_stand(lfeeg,lahee,lvb1,dsr);    % PHASE - burst first spikes
            aang_fs = [aang_fs paang_fs];
            ftm_fs0 = sum(exp(1).^(i*paang_fs)) / length(paang_fs);    % first trigonometric moment
            mvl_fs0 = abs(ftm_fs0);     % mean resultant length
            r_fs = [r_fs mvl_fs0];
            
            lsspo = sspo(sspo>ind1(k)&sspo<ind2(k)) - ind1(k);
            lsspo = round(lsspo/const);    % single spikes, downsample unit on 1000 Hz
            [paang_sp pdinx_sp] = laphase_stand(lfeeg,lahee,lsspo,dsr);    % PHASE - single spikes
            aang_sp = [aang_sp paang_sp];
            ftm_sp0 = sum(exp(1).^(i*paang_sp)) / length(paang_sp);    % first trigonometric moment
            mvl_sp0 = abs(ftm_sp0);     % mean resultant length
            r_sp = [r_sp mvl_sp0];
                        
            [paang_as pdinx_as] = laphase_stand(lfeeg,lahee,vdisc2,dsr);    % PHASE - all spikes
            aang_as = [aang_as paang_as];
            ftm_as0 = sum(exp(1).^(i*paang_as)) / length(paang_as);    % first trigonometric moment
            mvl_as0 = abs(ftm_as0);     % mean resultant length
            r_as = [r_as mvl_as0];
            
            lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
            lafsp = round(lafsp/const);    % all first spikes, downsample unit on 1000 Hz
            [paang_afsp pdinx_afsp] = laphase_stand(lfeeg,lahee,lafsp,dsr);    % PHASE - all first spikes
            aang_afsp = [aang_afsp paang_afsp];
            ftm_afsp0 = sum(exp(1).^(i*paang_afsp)) / length(paang_afsp);    % first trigonometric moment
            mvl_afsp0 = abs(ftm_afsp0);     % mean resultant length
            r_afsp = [r_afsp mvl_afsp0];
            
            libs = ibs(ibs>ind1(k)&ibs<ind2(k)) - ind1(k);
            libs = round(libs/const);    % intraburst spikes, downsample unit on 1000 Hz
            [paang_ibsang paang_sspoang paang_allang pcycnb] = laphaseb(lfeeg,lahee,lsspo,libs,lvb1,dsr);
            aang_ibsang = [aang_ibsang paang_ibsang];       % PHASE - "cycle first"
            aang_sspoang = [aang_sspoang paang_sspoang];
            aang_allang = [aang_allang paang_allang];
            cycnb = [cycnb pcycnb];
        end
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
end
if isempty(aang_fs)
    close all
    if exist('fname')
        disp([fname ' The frequency band is empty.'])
    else
        disp(['No file.'])
    end
    return
end

n_fs = length(aang_fs);     % burst first spikes
ftm_fs = sum(exp(1).^(i*aang_fs)) / n_fs;    % first trigonometric moment
ang_fs = angle(ftm_fs);   % mean angle
mvl_fs = abs(ftm_fs);     % mean resultant length
aang_fs = aang_fs * 180 / pi;
ang_fs = ang_fs * 180 / pi;
[nm_fs,xout_fs] = histc(aang_fs,edges);   % phase histogram
nm_fs = nm_fs(1:end-1);

n_sp = length(aang_sp);     % single spikes
ftm_sp = sum(exp(1).^(i*aang_sp)) / n_sp;    % first trigonometric moment
ang_sp = angle(ftm_sp);   % mean angle
mvl_sp = abs(ftm_sp);     % mean resultant length
aang_sp = aang_sp * 180 / pi;
ang_sp = ang_sp * 180 / pi;
[nm_sp,xout_sp] = histc(aang_sp,edges);   % phase histogram
nm_sp = nm_sp(1:end-1);

n_as = length(aang_as);     % all spikes
ftm_as = sum(exp(1).^(i*aang_as)) / n_as;    % first trigonometric moment
ang_as = angle(ftm_as);   % mean angle
mvl_as = abs(ftm_as);     % mean resultant length
aang_as = aang_as * 180 / pi;
ang_as = ang_as * 180 / pi;
[nm_as,xout_as] = histc(aang_as,edges);   % phase histogram
nm_as = nm_as(1:end-1);

n_afsp = length(aang_afsp);     % all first spikes
ftm_afsp = sum(exp(1).^(i*aang_afsp)) / n_afsp;    % first trigonometric moment
ang_afsp = angle(ftm_afsp);   % mean angle
mvl_afsp = abs(ftm_afsp);     % mean resultant length
aang_afsp = aang_afsp * 180 / pi;
ang_afsp = ang_afsp * 180 / pi;
[nm_afsp,xout_afsp] = histc(aang_afsp,edges);   % phase histogram
nm_afsp = nm_afsp(1:end-1);

n_ibsang = length(aang_ibsang);     % "cycle first"
ftm_ibsang = sum(exp(1).^(i*aang_ibsang)) / n_ibsang;    % first trigonometric moment
ang_ibsang = angle(ftm_ibsang);   % mean angle
mvl_ibsang = abs(ftm_ibsang);     % mean resultant length
aang_ibsang = aang_ibsang * 180 / pi;
ang_ibsang = ang_ibsang * 180 / pi;
[nm_ibsang,xout_ibsang] = histc(aang_ibsang,edges);   % phase histogram
nm_ibsang = nm_ibsang(1:end-1);
n_sspoang = length(aang_sspoang);
ftm_sspoang = sum(exp(1).^(i*aang_sspoang)) / n_sspoang;    % first trigonometric moment
ang_sspoang = angle(ftm_sspoang);   % mean angle
mvl_sspoang = abs(ftm_sspoang);     % mean resultant length
aang_sspoang = aang_sspoang * 180 / pi;
ang_sspoang = ang_sspoang * 180 / pi;
if ~isempty(aang_sspoang)
    [nm_sspoang,xout_sspoang] = histc(aang_sspoang,edges);   % phase histogram
    nm_sspoang = nm_sspoang(1:end-1);
else
    nm_sspoang =zeros(1,length(edges)-1);
    xout_sspoang = xout_ibsang;
end
n_allang = length(aang_allang);
ftm_allang = sum(exp(1).^(i*aang_allang)) / n_allang;    % first trigonometric moment
ang_allang = angle(ftm_allang);   % mean angle
mvl_allang = abs(ftm_allang);     % mean resultant length
aang_allang = aang_allang * 180 / pi;
ang_allang = ang_allang * 180 / pi;
[nm_allang,xout_allang] = histc(aang_allang,edges);   % phase histogram
nm_allang = nm_allang(1:end-1);
[nm_cycnb xout_cycnb] = hist(cycnb(cycnb>0),(1:10));     % distr. of burst no./cycle

figure(H1)
set(gcf,'Position',[62 378 1298 420])     % set figure window
subplot(1,2,1)      % all EPSPs
% bar(xout,nm_fs/length(aang_fs))
bar(cnts,nm_all'/n)
cmps = strread(fname,'%s','delimiter','_');
titlestr = [];
for tt = 1:length(cmps)
    titlestr = [titlestr ' ' cmps{tt}];
end
title(gca,[titlestr ' all EPSPs'])
x_lim = xlim;
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
% str = ['\it{Mean angle: }' '\bf ' num2str(ang)];
% text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
% text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
% str = ['\it{n: }' '\bf ' num2str(n)];
% text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

subplot(1,2,2)      % cycle first EPSPs
bar(cnts,nm_cfsp'/n_cfsp)
title(gca,[titlestr ' cycle first EPSPs'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
% str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
% text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
% text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
% str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
% text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

% Save
cd(resdir1)
fns = [fname(1:end-4) '_CYCLEFIRST.fig'];
saveas(H1,fns)
fns = [fname(1:end-4) '_CYCLEFIRST.jpg'];
saveas(H1,fns)
close all



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
vrs = version;
if isequal(vrs(1:5),'7.4.0') || isequal(vrs(1:5),'7.10.')
    files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
else
    files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
end
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function cyclen = eegfre(feeg,ahee,sr)

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 100 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
cl6 = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    lg = (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr);
    if ~lg
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
end
cyclen = mean(cl6) / sr * 1000;   % cycle length in ms

% -------------------------------------------------------------------------
function [ang inx] = laphase_stand(feeg,ahee,vdisc,sr)
%LAPHASE_STAND    Phase angles for unit relative to EEG.
%   [A I] = LAPHASE_STAND(FEEG,AHEE,VDISC,SR) calculates Hilbert phase 
%   angles (A) for discriminated unit (VDISC) relative to filtered EEG 
%   (FEEG), when sampling frequency is given in SR and Hilbert-transform of
%   the EEG in AHEE. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length.
%   Indices of discarded spikes of vdisc are returned in I.
%
%   See also HILBERT.

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 250 ms
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
inx = find(vdisc<fn(1));
plot(feeg)
plot(ahee,'k')
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
        plot(fn(k):fn(k+1),seeg,'r')
        plot(fn(k):fn(k+1),sahee,'r')
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);

% -------------------------------------------------------------------------
function [ang_ibsang ang_sspoang ang_allang cycnb] = laphaseb(feeg,ahee,sspo,ibs,vb1,sr)
%APHASEB    Phase angles for unit relative to EEG.
%   [IBSANG SSPOANG] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR) calculates
%   Hilbert phase angles for first intraburst and first single spike of
%   each cycle relative to the filtered EEG (FEEG), when sampling frequency
%   is given in SR, single spikes are given in SSPO, intraburst spikes are
%   given in IBS, burst first spikes are given in VB1 and Hilbert-transform
%   of the EEG in AHEE. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length. 
%
%   [IBSANG SSPOANG ALLANG] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR) returns
%   the phase of cycle first spikes as well.
%
%   [IBSANG SSPOANG ALLANG CYCNB] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR)
%   returns the number of bursts in each cycles.
%
%   See also HILBERT.

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 100 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
ibsang = [];
sspoang = [];
allang = [];
cycnb = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    if ~(axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        seeg = feeg(fn(k):fn(k+1));
        axs = max(seeg) - min(seeg);
        sahee = ahee(fn(k):fn(k+1));
        cycibs = ibs(ibs>fn(k)&ibs<fn(k+1));
        cycsspo = sspo(sspo>fn(k)&sspo<fn(k+1));
        cycvb1 = vb1(vb1>fn(k)&vb1<fn(k+1));
        if ~isempty(cycibs)
            if ~isempty(cycsspo)
                if cycibs(1) < cycsspo(1)
                    ibsang = [ibsang cycibs(1)];
                end
            else
                ibsang = [ibsang cycibs(1)];
            end
        end
        if ~isempty(cycsspo)
            if ~isempty(cycibs)
                if cycsspo(1) < cycibs(1)
                    sspoang = [sspoang cycsspo(1)];
                end
            else
                sspoang = [sspoang cycsspo(1)];
            end
        end
        if ~(isempty(cycibs) || isempty(cycsspo))
            allang = [allang min(cycibs(1),cycsspo(1))];
        elseif ~isempty(cycibs) && isempty(cycsspo)
            allang = [allang min(cycibs(1))];
        elseif isempty(cycibs) && ~isempty(cycsspo)
            allang = [allang min(cycsspo(1))];
        end
        cycnb(end+1) = length(cycvb1);
    end
end
ang_ibsang = ahee(ibsang);
ang_sspoang = ahee(sspoang);
ang_allang = ahee(allang);

% -------------------------------------------------------------------------
function [ang_allang cycnb] = aphaseb(eeg,vdisc,sr)
%APHASEB    Phase angles for unit relative to EEG.
%   [ALLANG CYCNB] = APHASEB(EEG,VDISC,SR) calculates Hilbert 
%   phase angles for cycle first events relative to the EEG, when sampling 
%   frequency is given in SR and event timestamps are given in VDISC. 
%   Cycles not fulfilling the following 2 criteria are discarded: (i) EEG
%   amp. higher then 2SD; (ii) min. 250 ms length. Number of events in each
%   cycle is returned in CYCNB.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 250 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
allang = [];
cycnb = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    if ~(axs < 2 * sd)  || (fn(k+1) - fn(k) < 0.25 * sr)
        seeg = feeg(fn(k):fn(k+1));
        axs = max(seeg) - min(seeg);
        sahee = ahee(fn(k):fn(k+1));
        cycvd = vdisc(vdisc>fn(k)&vdisc<fn(k+1));
        if ~isempty(cycvd)
            allang = [allang cycvd(1)];
        end
        cycnb(end+1) = length(cycvd);
    end
end
ang_allang = ahee(allang);