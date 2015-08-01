function aphase_corr
%APHASE_CORR   Correlation between EEG frequency and firing pattern.
%   APHASE_CORR calculates and saves correlation between 
%       (i) mean phase angle of all first spikes and EEG frequency
%       (ii) burstiness and EEG frequency.
%
%   See also APHASE and APHASERUN_BURST2.

% Directories
global DATAPATH
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl\';    % mat files
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Ketxyl\PhaseCorr\'];
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running APHASE_CORR...','Position',[360 250 275 50]);  %progress indicator
global WB
WB(end+1) = wb;

% Main
fre = zeros(1,sf);
meanang = zeros(1,sf);
brstness = zeros(1,sf);
for o = 1:sf
    fname = files_short{o};     % load
    ff = [inpdir1 fname];
    load(ff)
    eeg0 = data(:,2)';
    eeg = eeg0(1:20:end);    % downsample on 1000 Hz
    eeg2 = eeg0(1:100:end);    % downsample on 200 Hz
    clear data eeg0
    sr = 1000;
    ff2 = [inpdir2 fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    ssi = vdisc;       % allfirstspikes
    for k = 1:size(Burst,2)
        ssi(Burst(1,k)+1:Burst(2,k)) = 0;
    end
    afsp = ssi(ssi>0);
    afsp = round(afsp/20);    % downsample unit on 1000 Hz
    
    [aang_afsp dinx_afsp cyclen_afsp] = laphase(eeg,afsp,sr);    % PHASE - all first spikes
    n_afsp = length(aang_afsp);
    ftm_afsp = sum(exp(1).^(i*aang_afsp)) / n_afsp;    % first trigonometric moment
    ang_afsp = angle(ftm_afsp);   % mean angle
    mvl_afsp = abs(ftm_afsp);     % mean resultant length
    aang_afsp = aang_afsp * 180 / pi;
    ang_afsp = ang_afsp * 180 / pi;
    
    fre(o) = 1 / cyclen_afsp;
    meanang(o) = ang_afsp;
    brstness(o) = Burstiness;
    waitbar(o/sf)
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
text(0,2,num2str(Rangfre))
H2 = figure;
plot(brstness,fre,'.')
text(0.3,2,num2str(Rburstfre))

% Save
fns = [fname(1:end-4) '_ANGFRECORR.fig'];
saveas(H1,fns)
fns = [fname(1:end-4) '_BURSTFRECORR.fig'];
saveas(H1,fns)
save angfrecorr meanang fre Rangfre
save burstfrecorr brstness fre Rburstfre

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