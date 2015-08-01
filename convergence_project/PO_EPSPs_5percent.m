function [aang aang_cfsp resdir1] = PO_EPSPs_5percent(inpdir1)
%PO_EPSPS_5PERCENT   Phase analysis for EPSPs in PO neurons.
%   [A A_CFSP D] = PO_EPSPS_5PERCENT(DR) calculates phase histograms for
%   all EPSPs (A) and cycle first EPSPs (A_CFSP). Input directory (DR) is
%   specified in the caller function and results directory (D) is returned
%   to the same caller.
%
%   Note: slow oscillatoin cut at 500 ms (2 Hz); filter order is set to
%   2048 if data too short!
%
%   Call with AFRERESTRICTSTAND_CALL_LAYER5!
%
%   See also L5_SPIKING_5PERCENT and APHASERUN_BURST2_STAND_EPSP.

% Directories
global DATAPATH
% inpdir1 = [DATAPATH 'Hajni\EEGMPO\data\data for EEG_EPSP phase\mat\'];    % mat files
resdir1 = [DATAPATH 'Hajni_layer_5\PO_EPSPs_5percent\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir1);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running PO EPSPS 5PERCENT...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
[aang aang_cfsp] = main(inpdir1,resdir1,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function [aang aang_cfsp] = main(inpdir1,resdir1,files_short,sf,wb)

dsr = 1000;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
aang = [];
aang_cfsp = [];
for o = 1:sf   % segment cycle
    fname = files_short{o};     % load
    ff = [inpdir1 fname];
    load(ff)
    sr = 1 / EEG1.interval;
    const = sr / dsr;
    eeg0 = EEG1.values';
    eeg = eeg0(1:const:end);    % downsample on 1000 Hz
    vdisc = EPSPs.times';
    vdisc = round(vdisc*dsr);    % sample on 1000 Hz
    clear EEG1 EPSPs eeg0
    
    [paang dinx] = laphase_stand(eeg,vdisc,dsr);    % PHASE - all EPSPs
    aang = [aang paang];
    [paang_cfsp cycnb] = aphaseb(eeg,vdisc,dsr);    % PHASE - cycle first EPSPs
    aang_cfsp = [aang_cfsp paang_cfsp];
    waitbar(o/sf,wb);
end
    
% Phase histograms
n = length(aang);
ftm = sum(exp(1).^(i*aang)) / n;    % first trigonometric moment
ang = angle(ftm);   % mean angle
mvl = abs(ftm);     % mean resultant length
aang = aang * 180 / pi;
ang = ang * 180 / pi;
[nm_all,xout_all] = histc(aang,edges);   % phase histogram
nm_all = nm_all(1:end-1);

n_cfsp = length(aang_cfsp);
ftm_cfsp = sum(exp(1).^(i*aang_cfsp)) / n_cfsp;    % first trigonometric moment
ang_cfsp = angle(ftm_cfsp);   % mean angle
mvl_cfsp = abs(ftm_cfsp);     % mean resultant length
aang_cfsp = aang_cfsp * 180 / pi;
ang_cfsp = ang_cfsp * 180 / pi;
[nm_cfsp,xout_cfsp] = histc(aang_cfsp,edges);   % phase histogram
nm_cfsp = nm_cfsp(1:end-1);

% Display percentiles
pct5 = prctile(aang_cfsp,5);
disp(['5%: ' num2str(pct5)])
pct10 = prctile(aang_cfsp,10);
disp(['10%: ' num2str(pct10)])

% Plot
H1 = figure;
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
str = ['\it{Mean angle: }' '\bf ' num2str(ang)];
text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n)];
text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

subplot(1,2,2)      % cycle first EPSPs
bar(cnts,nm_cfsp'/n_cfsp)
title(gca,[titlestr ' cycle first EPSPs'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

H2 = figure;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm_cfsp,xout_cfsp] = histc(aang_cfsp,edges);   % phase histogram
nm_cfsp = nm_cfsp(1:end-1);plot(cnts,nm_cfsp'/n_cfsp,'k')
hold on
ml5 = prctile(aang_cfsp,10);
xml5 = linterp(cnts,nm_cfsp'/n_cfsp,ml5);
area([cnts(cnts<ml5) ml5],[nm_cfsp(cnts<ml5)/n_cfsp xml5],'FaceColor','red')
title(gca,'Cycle first EPSPs')
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

H3 = figure;
edges = -180:10:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm_cfsp,xout_cfsp] = histc(aang_cfsp,edges);   % phase histogram
nm_cfsp = nm_cfsp(1:end-1);plot(cnts,nm_cfsp'/n_cfsp,'k')
hold on
ml5 = prctile(aang_cfsp,10);
xml5 = linterp(cnts,nm_cfsp'/n_cfsp,ml5);
area([cnts(cnts<ml5) ml5],[nm_cfsp(cnts<ml5)/n_cfsp xml5],'FaceColor','red')
title(gca,'Cycle first EPSPs')
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

% Save
cd(resdir1)
fns = [fname(1:end-4) '_CYCLEFIRST.fig'];
saveas(H1,fns)
fns = [fname(1:end-4) '_CYCLEFIRST.jpg'];
saveas(H1,fns)
fns = [fname(1:end-4) '_CYCLEFIRSTb.fig'];
saveas(H2,fns)
fns = [fname(1:end-4) '_CYCLEFIRSTb.jpg'];
saveas(H2,fns)
fns = [fname(1:end-4) '_CYCLEFIRSTc.fig'];
saveas(H3,fns)
fns = [fname(1:end-4) '_CYCLEFIRSTc.jpg'];
saveas(H3,fns)
save([fname(1:end-4) 'CYCLEFIRST.mat'],'aang','aang_cfsp','pct5','pct10')
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
function [ang inx] = laphase_stand(eeg,vdisc,sr)
%APHASE_STAND    Phase angles for unit relative to EEG.
%   [A I] = LAPHASE_STAND(EEG,VDISC,SR) calculates Hilbert phase angles (A)
%   for discriminated unit (VDISC) relative to standardized EEG, when
%   sampling frequency is given in SR. Cycles not fulfilling the following
%   2 criteria are discarded: (i) EEG amp. higher then 2SD; (ii) min. 500
%   ms length. Indices of disclosed spikes of vdisc are returned in I.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
try
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,eeg);
catch
    flt = fir1(2048,5/nqf,'low');
    disp('Filter order: 2048')
    feeg = filtfilt(flt,1,eeg);
end
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 500 ms
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
inx = find(vdisc<fn(1));
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.5 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);



% -------------------------------------------------------------------------
function [ang_allang cycnb] = aphaseb(eeg,vdisc,sr)
%APHASEB    Phase angles for unit relative to EEG.
%   [ALLANG CYCNB] = APHASEB(EEG,VDISC,SR) calculates Hilbert 
%   phase angles for cycle first events relative to the EEG, when sampling 
%   frequency is given in SR and event timestamps are given in VDISC. 
%   Cycles not fulfilling the following 2 criteria are discarded: (i) EEG
%   amp. higher then 2SD; (ii) min. 500 ms length. Number of events in each
%   cycle is returned in CYCNB.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
try
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,eeg);
catch
    flt = fir1(2048,5/nqf,'low');
    disp('Filter order: 2048')
    feeg = filtfilt(flt,1,eeg);
end
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 500 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
allang = [];
cycnb = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    if ~(axs < 2 * sd)  || (fn(k+1) - fn(k) < 0.5 * sr)
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