function [aang aang_cfsp resdir1] = L5_spiking_5percent(inpdir1)
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
wb = waitbar(0,'Running L5 SPIKING 5PERCENT...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
[aang aang_cfsp] = main(inpdir1,resdir1,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function [aang aang_cfsp] = main(inpdir1,resdir1,files_short,sf,wb)

sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
aang = [];
aang_cfsp = [];
seglen = 0;
% H1 = figure;
mlen = 94;  % calculate results from 94 seconds
for o = 1:sf
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
        
    [paang dinx] = aphase_stand(eeg,vdisc,dsr);    % PHASE - all EPSPs
    aang = [aang paang];
    [paang_cfsp cycnb] = aphaseb(eeg,vdisc,dsr);    % PHASE - cycle first EPSPs
    aang_cfsp = [aang_cfsp paang_cfsp];
    waitbar(o/sf,wb);
    if seglen > mlen
        break
    end
end
% n = length(aang);
% ftm = sum(exp(1).^(i*aang)) / n;    % first trigonometric moment
% ang = angle(ftm);   % mean angle
% mvl = abs(ftm);     % mean resultant length
% aang = aang * 180 / pi;
% ang = ang * 180 / pi;
% [nm_all,xout_all] = histc(aang,edges);   % phase histogram
% nm_all = nm_all(1:end-1);
% 
% n_cfsp = length(aang_cfsp);
% ftm_cfsp = sum(exp(1).^(i*aang_cfsp)) / n_cfsp;    % first trigonometric moment
% ang_cfsp = angle(ftm_cfsp);   % mean angle
% mvl_cfsp = abs(ftm_cfsp);     % mean resultant length
% aang_cfsp = aang_cfsp * 180 / pi;
% ang_cfsp = ang_cfsp * 180 / pi;
% [nm_cfsp,xout_cfsp] = histc(aang_cfsp,edges);   % phase histogram
% nm_cfsp = nm_cfsp(1:end-1);
% 
% figure(H1)
% set(gcf,'Position',[62 378 1298 420])     % set figure window
% subplot(1,2,1)      % all spikes
% % bar(xout,nm_fs/length(aang_fs))
% bar(cnts,nm_all'/n)
% cmps = strread(fname,'%s','delimiter','_');
% titlestr = [];
% for tt = 1:length(cmps)
%     titlestr = [titlestr ' ' cmps{tt}];
% end
% title(gca,[titlestr ' all spikes'])
% x_lim = xlim;
% y_lim = ylim;
% axis([-200 200 y_lim(1) y_lim(2)])
% % str = ['\it{Mean angle: }' '\bf ' num2str(ang)];
% % text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% % str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
% % text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
% % str = ['\it{n: }' '\bf ' num2str(n)];
% % text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% set(gca,'XTick',[-180 -90 0 90 180])
% 
% subplot(1,2,2)      % cycle first spikes
% bar(cnts,nm_cfsp'/n_cfsp)
% title(gca,[titlestr ' cycle first spikes'])
% y_lim = ylim;
% axis([-200 200 y_lim(1) y_lim(2)])
% % str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
% % text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% % str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
% % text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
% % str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
% % text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% set(gca,'XTick',[-180 -90 0 90 180])
% 
% % Save
% cd(resdir1)
% fns = [fname(1:end-4) '_CYCLEFIRST.fig'];
% saveas(H1,fns)
% fns = [fname(1:end-4) '_CYCLEFIRST.jpg'];
% saveas(H1,fns)
% close all



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