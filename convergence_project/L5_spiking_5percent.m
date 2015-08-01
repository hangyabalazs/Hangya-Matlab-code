function [aang aang_cfsp resdir1] = L5_spiking_5percent(inpdir1)
%L5_SPIKING_5PERCENT   Phase analysis for layer 5 pyramidal neurons.
%   [A A_CFSP D] = L5_SPIKING_5PERCENT(DR) calculates phase histograms for
%   all spikes (A) and cycle first spikes (A_CFSP). Input directory (DR) is
%   specified in the caller function and results directory (D) is returned
%   to the same caller. Segment are cut at a predefined time instance to
%   provide constituents for a time-balanced pooled sample.
%
%   Note: slow oscillatoin cut at 500 ms (2 Hz); filter order is set to
%   2048 if data too short!
%
%   See also AFRE_RESTRICT_STAND_LAYER5.

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
wb = waitbar(0,'Running L5 SPIKING 5PERCENT...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
[aang aang_cfsp] = main(inpdir1,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function [aang aang_cfsp] = main(inpdir1,files_short,sf,wb)

sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
aang = [];
aang_cfsp = [];
seglen = 0;
% H1 = figure;
mlen = 137;  % calculate results from 94 seconds
for o = 1:sf   % segment cycle
    fname = files_short{o}     % filename
    cmps = strread(fname(1:end-4),'%s','delimiter','_');
    ff = [inpdir1 fname];       % load
    load(ff)
    seglen0 = seglen;
    lseglen = str2double(cmps{5}) - str2double(cmps{4});
    seglen = seglen + lseglen;
    if seglen > mlen
        restrictseg = mlen - seglen0;
        seglen0 + restrictseg
        eeg = eeg(1:restrictseg*sr);
    end
    eeg = eeg(1:const:end);    % sample at 1000 Hz
    vdisc = round(vdisc/const);
        
    [paang dinx] = laphase_stand(eeg,vdisc,dsr);    % PHASE - all EPSPs
    aang = [aang paang];
    [paang_cfsp cycnb] = aphaseb(eeg,vdisc,dsr);    % PHASE - cycle first EPSPs
    aang_cfsp = [aang_cfsp paang_cfsp];
    waitbar(o/sf,wb);
    if seglen > mlen
        break
    end
end



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