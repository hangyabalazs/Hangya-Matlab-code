function afre_restrict_bic_stand(inpdir1)
%AFRE_RESTRICT_BIC_STAND   Phase and burst analysis restricted to a frequency band.
%   AFRE_RESTRICT_BIC_STAND(DR) calculates and saves phase and burst
%   analysis results for a given EEG frequency band (1.1 - 1.6 Hz) for
%   baseline-bicuculline comparison. See ACLUSTERCUT and APHASE_BURST2_STAND
%   for details on the analysis. Input directory should be given as an
%   argument (DR). EEG is standardized for phase calculations.
%
%   See also ACLUSTERCUT and APHASERUN_BURST2_STAND.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
inpdir_bas = [inpdir1 'bas\'];
inpdir_bic = [inpdir1 'bic\'];
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_burst_bic_stand\'];
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
[wb,awb1,awb2] = waitbar2([0 0],'Running AFRE RESTRICT BIC STAND...');
global WB
WB(end+1) = wb;

% Main
main(inpdir1,inpdir2,resdir1,files_short_bas,sf_bas,wb,'bas');
main(inpdir1,inpdir2,resdir1,files_short_bic,sf_bic,wb,'bic');

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir1,inpdir2,resdir1,files_short,sf,wb,bob);

sr = 20000;
dsr = 1000;
const = sr / dsr;
preburst = [];
burstiness = [];
intraburstfreq = [];
ibspno = [];
burstlength = [];
burstfreq = [];
frate = [];
isimtx = {};
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
    
    nqf = dsr / 2;      % filtering EEG
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,eeg(1:const:end));
    feeg = (feeg - mean(feeg)) / std(feeg);
    ahee = angle(hilbert(feeg));    % Hilbert-transformation
    
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
        
% Restrict frequency band
        eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/const);
        
        cyclen = eegfre(lfeeg,lahee,dsr);    % filter EEG, Hilbert-transform
        freq = 1 / cyclen * 1000;
        if 1.1 < freq & freq < 1.6        % restrict frequency band
            
% Burst statistics
            [lvb lburst] = locvburst(vburst,Burst,ind1(k),ind2(k));
            preburst = [preburst lburst];
            burstnum = size(lvb,2);
            intraburstiv = [];
            intraburstnum = zeros(1,burstnum);
            for j = 1:burstnum      % computing intraburstiv
                b = vdisc(vdisc>=lvb(1,j)&vdisc<=lvb(2,j));
                db = diff(b);
                intraburstiv = [intraburstiv db];
                intraburstnum(j) = length(b);   %intraburst spike number
            end
            burstiness = [burstiness (length(intraburstiv) + burstnum) / length(vd)];
            burstlen = (lvb(2,:) - lvb(1,:)) / sr;
            burstlength = [burstlength burstlen];
            if ~isempty(intraburstnum)
                intraburstfreq = [intraburstfreq (intraburstnum-1)./burstlen];
            end
            ibspno = [ibspno intraburstnum];
            burstfreq = [burstfreq sr * (burstnum - 1)  / (lvb(2,end) - lvb(1,1))];
            efflen = (vd(end) - vd(1)) / sr;
            frate = [frate (length(vd) - 1) / efflen];
        end
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
end
if isempty(frate)
    close all
    if exist('fname')
        disp([fname ' The frequency band is empty.'])
    else
        disp(['No ' bob 'file.'])
    end
    return
end

Burstiness = burstiness;     % burst parameters
IntraBurstFrequency.mean = mean(intraburstfreq);
IntraBurstFrequency.sd = std(intraburstfreq);
IntraBurstFrequency.all = intraburstfreq;
IntraBurstSpikeNumber.mean = mean(ibspno);
IntraBurstSpikeNumber.sd = std(ibspno);
IntraBurstSpikeNumber.all = ibspno;
BurstLength.mean = mean(burstlength);
BurstLength.sd = std(burstlength);
BurstLength.all = burstlength;
BurstFrequency = burstfreq;
FiringRate = frate;
Burst = preburst;

% Save
cd(resdir1)
fn = [fname(1:end-4) '_CLUST2.mat'];
save(fn,'Burstiness','IntraBurstFrequency','IntraBurstSpikeNumber','BurstLength',...
    'BurstFrequency','FiringRate','Burst','vdisc')



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
function [lvburst lburst] = locvburst(vburst,Burst,ind1,ind2)

fst = vburst(1,:);
lst = vburst(2,:);
fi = find(fst>ind1,1,'first');
li = find(lst<ind2,1,'last');
lvburst = vburst(:,fi:li);
lburst = Burst(:,fi:li);

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