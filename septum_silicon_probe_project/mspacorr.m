function mspacorr
%MSPACORR   Autocorrelation.
%   MSPACORR imports unit and EEG data of the MSHCsp project; it calculates
%   unit autocorrelations for theta and non-theta segments.
%
%   See also CZACORR, RAPHEVIEW2 and AFRE_RESTRICT_SILICON.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'MSHCsp\Viktor5\'];
resdir = [DATAPATH 'MSHCsp\Acorr\Viktor5\'];
thetadir = [DATAPATH 'MSHCsp\Wavelet\theta_segments\'];
nonthetadir = [DATAPATH 'MSHCsp\Wavelet\nontheta_segments\'];
mm = pwd;

% Filelist
% rid = '1234';
% flist = {'201003171' '201003172' '201003173' '201003174'};
rid = '1237';
flist = {'201003291' '201003292' '201003293' '201003294'};
sf = length(flist);

% Main
sr = 20000;
% dsr = 1000;
% const = sr / dsr;
% edges = -180:20:180;     % edges for phase histogram
% cnts = (edges(1:end-1) + edges(2:end)) / 2;
for o = 1:sf
    fname = flist{o}     % 'filename'
%     ff = [inpdir 'EEG\' rid '_' fname '_eeg.mat'];       % load EEG
%     load(ff)
%     len = length(eeg);
    
    ff = [thetadir 'THETA_SEGMENTS_' rid '_' fname '.mat'];     % load theta segments
    load(ff)
    ThetaSegments = uniteseg(ThetaSegments,sr);     % drop gaps < 0.5 s
    ThetaSegments = short_killer(ThetaSegments);    % drop segments < 3 s
    lent = ThetaSegments(2,:) - ThetaSegments(1,:);
    minx = find(lent==max(lent));
    th1 = ThetaSegments(1,minx);
    th2 = ThetaSegments(2,minx);
    
    ff = [nonthetadir 'NONTHETA_SEGMENTS_' rid '_' fname '.mat'];     % load non-theta segments
    load(ff)
    NonThetaSegments = short_killer(NonThetaSegments);    % drop segments < 3 s
    lenn = NonThetaSegments(2,:) - NonThetaSegments(1,:);
    minx = find(lenn==max(lenn));
    nth1 = NonThetaSegments(1,minx);
    nth2 = NonThetaSegments(2,minx);
    
    for shankno = 1:4
        ff = [inpdir fname '.res.' num2str(shankno)];       % load unit
        res = textread(ff);     % spike timestamps
        ff = [inpdir fname '.clu.' num2str(shankno)];
        if ~b_isfilename(ff)
            continue
        end
        clu = textread(ff);     % indices of sorted cells
        clu = clu(2:end);
        ncells = max(clu);        % number of sorted cells (including multiunit)
        
%         dind = round(1:const:len);      % downsample EEG on 1000 Hz
%         deeg = eeg(dind);
%         nqf = dsr / 2;      % filtering EEG
%         flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
%         feeg = filtfilt(flt,1,deeg);
%         ahee = angle(hilbert(feeg));    % Hilbert-transformation

        for nc = 2:ncells
            vdisc = res(clu==nc)';
            vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
            if ~isempty(vdisc_theta)
                ac = racorr(vdisc_theta/sr,500);    % theta autocorrelation
                H = figure;
                bar(linspace(-1000,1000,length(ac)),ac)
                fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_THACORR.fig'];
                saveas(H,fns)
            end
            vdisc_nontheta = vdisc(vdisc>nth1&vdisc<nth2) - nth1;
            if ~isempty(vdisc_nontheta)
                ac = racorr(vdisc_nontheta/sr,500);    % non-theta autocorrelation
                H = figure;
                bar(linspace(-1000,1000,length(ac)),ac)
                fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_NTHACORR.fig'];
                saveas(H,fns)
            end
        end
    end
end
cd(mm)

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
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);



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