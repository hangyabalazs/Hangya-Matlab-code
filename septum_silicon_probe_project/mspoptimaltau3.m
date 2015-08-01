function mspoptimaltau3(ratno)
%MSPOPTIMALTAU3   TE calculation.
%   MSPOPTIMALTAU3 calculates transfers entropy for the all theta segments
%   either spontane or induced and longest non-theta. It uses a predefined
%   value for optimal tau.
%
%   See also MSPTRTE, MSPTHETASELECTORRUN5.

%   matlabpool local 7

% Input argument check
% error(nargchk(1,1,nargin))
ratno = 8;

% Directories
global DATAPATH
global DATADIR
% for ratno = 7:9
no = num2str(ratno);
inpdir = [DATADIR 'MSHCsp\viktor' no '\'];
segdir = [DATAPATH 'MSHCsp\analysis\segments\viktor' no '\'];
try
    tedir = [DATAPATH 'MSHCsp\analysis\opttau\viktor' no '\'];
    cd(tedir)
catch
    mkdir(tedir)
end
mm = pwd;

% Filelist
load([DATADIR 'MSHCsp\flist.mat' ])
flist = flist{ratno};
sf = length(flist);

% Optimal tau
% load or type optima tau here
% name it 'tau'

% Main
sr = 20000;    % sampling rate in ms
dsr = 400;
const = sr / dsr;
mm = ['s'; 'p'; 'n'];
for (o =  1:sf)    % no. of the measurement
    fname = num2str(flist{o});
    disp(fname)
    for ch = 1:32 % mspreadeeg reads 33th ch as 1st
        disp(ch)
        
        % Load theta/nontheta segments
        tsg = load([segdir  'segments_' fname '_' num2str(ch) '.mat']);
        segdat = tsg;
        disp([segdir  'segments_' fname '_' num2str(ch) '.mat']);
        for tpp =  1:3             % spontane (1), pinch-induced (2) theta, non-theta(3)
            segtype = mm(tpp);
            disp(segtype)
            if tpp == 3
                segtype = 'Non';
            end
            segname = [segtype 'ThetaSegments'];
%             segs = ceil(eval([segtype 'ThetaSegments']));
            segs = getfield(segdat,segname);
            for segno = 1:size(segs,2) % analyze each segments
                disp(['segment number: ' num2str(segno)])
                if tpp == 3 && segno > 5 % 5 non-theta segments are enough
                    break
                end
                seg = segs(1,segno):segs(2,segno);
                seglen = length(seg);
                
                % Load EEG
                ffnm = [inpdir 'EEG\' fname '_eeg.mat'];
                eeg = [];
                teeg = load(ffnm);
                eeg = teeg.eeg;
                disp(ffnm)
                eegall= eeg;
                eeg = eegall{ch};
                disp(['eeg length: ', num2str(length(eeg))])
                disp([num2str(seg(1)) '_' num2str(seg(end))])
                if seg(end) > length(eeg)
                    srpl = seg(end) - length(eeg);
                    seg = seg(1:seglen-srpl);
                end
                eeg_seg = (eeg(seg))';
                disp('If this appears, a critical point has benn passed')
                
                % Load unit
                for shankno = 1:4
                    ff = [inpdir fname '.res.' num2str(shankno)];       % load unit
                    res = textread(ff);     % spike timestamps
                    ff = [inpdir fname '.clu.' num2str(shankno)];
                    clu = textread(ff);     % indices of sorted cells
                    clu = clu(2:end);
                    ncells = max(clu);        % number of sorted cells (including multiunit)
                    
                    for nc = 1:ncells
                        vdisc = res(clu==nc)';
                        vdisc = round(vdisc/const);  % resample at 'dsr'
                        vdisc_seg = vdisc(ismember(vdisc,seg)==1) - seg(1) + 1;
                        if length(vdisc_seg) < 20    % segments with min. 20 spike
                            break
                        end
                        
                        % TRANSFER ENTROPY
                        % Entropy calculation
                        
                        % Wavelet transformation of the EEG
                        abs1 = eeg_wavelet(eeg_seg,dsr);
                        
                        % Create random unit (shuffle ISIs)
                        psvd_seg = isi_shuffle(vdisc_seg);
                        
                        % Sinc convolution & wavelet transformation of unit
                        abs2 = unit_wavelet(vdisc_seg,seglen,sr,dsr);   % unit
                        abs3 = unit_wavelet(psvd_seg,seglen,sr,dsr);    % random unit
                        
                        if ~isequal(size(abs1,2),size(abs2,2))    % adjust size
                            mininx = min([size(abs1,2) size(abs2,2)]);
                            abs1 = abs1(:,1:mininx);
                            abs2 = abs2(:,1:mininx);
                            abs3 = abs3(:,1:mininx);
                        end
                        
                        min_abs1 = min(abs1(:));    % calculate extrema
                        max_abs1 = max(abs1(:));
                        min_abs2 = min(abs2(:));
                        max_abs2 = max(abs2(:));
                        min_abs3 = min(abs3(:));
                        max_abs3 = max(abs3(:));
                        
                        cTE = struct('seg',struct('eu',[],'ue',[]));
                        cNTE = struct('seg',struct('eu',[],'ue',[]));
                        cBias = struct('seg',struct('eu',[],'eu_shuffled',[],'ue',[],'ue_shuffled', []));
                        [TE TE_corr TE_bias H_X2FcX2] = ltren(abs1,abs2,...
                            min_abs1,max_abs1,min_abs2,max_abs2,dsr,tau);  % theta, EEG->unit
                        TE_seg_eu = TE;
                        TE_seg_eu_bias = TE_bias;
                        [TE TE_corr TE_bias] = ltren(abs1,abs3,...
                            min_abs1,max_abs1,min_abs3,max_abs3,dsr,tau);
                        TE_seg_eu_shuffled = TE;
                        TE_seg_eu_shuffled_bias = TE_bias;
                        NTE_seg_eu = (TE_seg_eu - TE_seg_eu_shuffled) / H_X2FcX2;
                        NTE_seg_eu = max(NTE_seg_eu,0);
                        
                        [TE TE_corr TE_bias H_X2FcX2] = ltren(abs2,abs1,...
                            min_abs2,max_abs2,min_abs1,max_abs1,dsr,tau);  % theta, unit->EEG
                        TE_seg_ue = TE;
                        TE_seg_ue_bias = TE_bias;
                        [TE TE_corr TE_bias] = ltren(abs3,abs1,...
                            min_abs3,max_abs3,min_abs1,max_abs1,dsr,tau);
                        TE_seg_ue_shuffled = TE;
                        TE_seg_ue_shuffled_bias = TE_bias;
                        NTE_seg_ue = (TE_seg_ue - TE_seg_ue_shuffled) / H_X2FcX2;
                        NTE_seg_ue = max(NTE_seg_ue,0);
                        
                        if abs(NTE_seg_ue) + abs(NTE_seg_eu) == Inf
                            keyboard
                        end
                        
                        DF_seg_eu = (NTE_seg_eu - NTE_seg_ue) / (NTE_seg_eu + NTE_seg_ue);    % direction of flow (pos: EEG->unit)
                        
                        cTE.seg.eu = TE_seg_eu;
                        cTE.seg.ue = TE_seg_ue;
                        cNTE.seg.eu = NTE_seg_eu;
                        cNTE.seg.ue = NTE_seg_ue;
                        cBias.seg.eu = TE_seg_eu_bias;
                        cBias.seg.eu_shuffled = TE_seg_eu_shuffled_bias;
                        cBias.seg.ue = TE_seg_ue_bias;
                        cBias.seg.ue_shuffled = TE_seg_ue_shuffled_bias;
                        cDF_seg_eu = DF_seg_eu;
                        cNTE_seg_eu = NTE_seg_eu;
                        cNTE_seg_ue = NTE_seg_ue;
                            
                        mdf_seg = cDF_seg_eu;
                        mDF_seg(o) = mdf_seg;
                        
                        % Save
                        mfn = [tedir fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc)...
                            '_' mm(tpp) num2str(segno) '_mTE.mat' ];
                        parlsave(mfn,'mdf_seg')
                        close all
                        
                        % TE significance values
                        snm = 100;
                        sTE_seg_eu = zeros(1,snm);
                        sTE_seg_ue = zeros(1,snm);
                        for k = 1:snm
                            psvd_seg  = isi_shuffle(vdisc_seg); % create random unit (shuffle ISIs)
                            
                            abs4 = unit_wavelet(psvd_seg, seglen,sr,dsr); % sinc convolution & wavelet transformation of random unit
                            if ~isequal(size(abs1,2),size(abs4,2))    % adjust size
                                mininx2 = min([size(abs1,2) size(abs4,2)]);
                                abs1 = abs1(:,1:mininx2);
                                abs4 = abs4(:,1:mininx2);
                            end
                            min_abs4 = min(abs4(:));
                            max_abs4 = max(abs4(:));
                            
                            TE = ltren2(abs1,abs4,...
                                min_abs1,max_abs1,min_abs4,max_abs4,dsr,mtau_seg);  % theta, EEG->unit
                            sTE_seg_eu(k) = TE;
                            
                            TE = ltren2(abs4,abs1,...
                                min_abs4,max_abs4,min_abs1,max_abs1,dsr,mtau_seg);  % theta, unit->EEG
                            sTE_seg_ue(k) = TE;
                        end
                        mfn = [tedir fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc) '_' mm(tpp) num2str(segno) '_TEsign.mat' ];
                        parlsave2(mfn,'sTE_seg_eu','sTE_seg_ue');
                    end
                end
            end
        end
    end
end
cd(mm)

% -------------------------------------------------------------------------
% RANDOM UNIT
% -------------------------------------------------------------------------
function [psvd_seg] = isi_shuffle(vdisc_seg)

% Interspike intervals
isi_seg = diff(vdisc_seg);


% Randomize segment
lit = length(isi_seg);    % shuffle ISIs
rp = randperm(lit);
while any(rp==(1:lit))
    rp = randperm(lit);
end
psi1 = [];
for it = 1:lit
    psi1 = [psi1 isi_seg(rp(it))];
end
psvd_seg = [vdisc_seg(1) vdisc_seg(1)+cumsum(psi1)];




% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function pow = eeg_wavelet(dat,dsr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / dsr;
pad = 1;
dj = 0.08;
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt;
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;

fnd = find(f>6);    % frequency band bounderies
pwind1 = fnd(end);
fnd = find(f<2.5);
pwind2 = fnd(1);
f = f(pwind1:pwind2);
fperiod = 1 ./ f;
s = fperiod ./ c;

param = -1;
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = wavelet(dat,dt,pad,s,mother,param);
pow = abs(wave).^2;



% -------------------------------------------------------------------------
function pow = unit_wavelet(vdisc,seglen,sr,dsr)

% Sinc convolution
vdisco = vdisc * (sr / dsr);
fs = sr;     % unit
dto = 1 / fs;
fcut = 100;
fsnew = dsr;
dtnew = 1 / fsnew;
fsratio = fsnew / fs;
told = vdisco * dto * fcut;
lenu = seglen * (sr / fsnew);
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);

zint = 0;
for k = 1:lentold
    %                             zint = zint + sinc(tnew - told(k));
    zint = zint + sinc(tnew - told(k));
end

%                         tnew2 = (-lenu*fsratio:lenu*fsratio) * dtnew * fcut;
%                         psinc = sinc(tnew2);
%                         lt = length(tnew);
%                         told2 = round(told/dtnew/fcut);
%                         zint = zeros(1,lt);
%                         for k = 1:lentold
%                             zint = zint + psinc(lt-told2(k)+2:2*lt-told2(k)+1);
%                         end

% Prepare for wavelet transformation
variance = std(zint) ^ 2;
zint = (zint - mean(zint)) / sqrt(variance) ;
n = length(zint);
dt = 1 / dsr;
pad = 1;
dj = 0.08;
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt;
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;

fnd = find(f>6);    % frequency band bounderies
pwind1 = fnd(end);
fnd = find(f<2.5);
pwind2 = fnd(1);
f = f(pwind1:pwind2);
fperiod = 1 ./ f;
s = fperiod ./ c;

param = -1;
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = wavelet(zint,dt,pad,s,mother,param);
pow = abs(wave).^2;



% -------------------------------------------------------------------------
function [wave,period,scale,coi] = wavelet(Y,dt,pad,s,mother,param)
%WAVELET   Wavelet with scales for every frequency.
%
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

n1 = length(Y);

% Construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
    base2 = fix(log(n1)/log(2)+0.4999);   % power of 2 nearest to N
    x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

% Construct wavenumber array used in transform [Eqn(5)]
k = 1:fix(n/2);
k = k .* ((2 .* pi) / (n * dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

% Compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

% Construct SCALE array & empty PERIOD & WAVE arrays
scale = s;

% Loop through all scales and compute transform
wave = zeros(length(scale),n1);  % define the wavelet array
for a1 = 1:length(scale)
    [daughter,fourier_factor,coi,dofmin] = wave_bases2(mother,k,scale(a1),param);
    ifd = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
    ifd = ifd(1:n1);
    wave(a1,:) = ifd;
end

period = fourier_factor * scale;
coi = coi * dt * [1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]



% -------------------------------------------------------------------------
% ENTROPY
% -------------------------------------------------------------------------
function [TE TE_corr TE_bias H_X2FcX2] = ltren(W1,W2,minW1,maxW1,minW2,maxW2,sr,tau)
%TREN   Transfer entropy.
%   TE = TREN(W1,W2,SR) calculates transfer entropy for time series W1 and
%   W2 sampled on SR (TE_W1->W2).
%
%   TE = TREN(W1,W2,SR,TAU) uses TAU ms as time lag between current and
%   "future" values.
%
%   [TE TE_CORR TE_BIAS] = TREN(W1,W2,SR,TAU) returns unbiased estimate of
%   transfer entropy applying Treves-Panzeri bias correction algorithm. The
%   amount of the bias of the original estimate is also returned.
%
%   [TE TE_CORR TE_BIAS H] = TREN(W1,W2,SR,TAU) returns H(X2F|X2)
%   conditonal entropy for transfer entropy normalization.
%
%   References:
%   (1) Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007)
%   Correcting for the sampling bias problem in spike train information
%   measures. J Neurophysiol 98:1064-1072.
%   (2) Gour�vitch B, Eggermont JJ (2007) Evaluating information transfer
%   between auditory cortical neurons. J Neurophysiol 97:2533-2543.
%   (3) Imas OA, Ropella KM, Ward BD, Wood JD, Hudetz AG (2005) Volatile
%   anesthetics disrupt frontal-posterior recurrent information transfer at
%   gamma frequencies in rat. Neurosci Lett 387:145-50.

% Input argument check
error(nargchk(7, 8, nargin))
if nargin < 8
    tau = 100;   % time lag in ms
end

% Extracting the random variables
tau2 = tau / 1000 * sr;     % time lag in data points
X1 = W1(:,1:end-tau2);
X2 = W2(:,1:end-tau2);
X2F = W2(:,tau2+1:end);
n = numel(X1);

% Calculating joint histograms
bno = fix(exp(0.626+0.4*log(n-1)));   % bin number for histogram estimates
bno = 30;
minX1 = minW1;    % discretization of X1
maxX1 = maxW1;
binwidth1 = (maxX1 - minX1) ./ bno;
xx1 = minX1 + binwidth1 * (0:bno);   % bin edges
xx1(length(xx1)) = maxX1;
xx1(1) = -inf;
nbin1 = length(xx1);

minX2 = minW2;    % discretization of X2
maxX2 = maxW2;
binwidth2 = (maxX2 - minX2) ./ bno;
xx2 = minX2 + binwidth2 * (0:bno);   % bin edges
xx2(length(xx2)) = maxX2;
xx2(1) = -inf;
nbin2 = length(xx2);

h_X2F_X2_X1 = zeros(nbin2-1,nbin2-1,nbin1-1);   % P(X2F, X2, X1)
t1 = X2F(:) - minX2;
t2 = X2(:) - minX2;
t3 = X1(:) - minX1;
p1 = fix((t1-1000000*eps)/binwidth2) + 1;
p2 = fix((t2-1000000*eps)/binwidth2) + 1;
p3 = fix((t3-1000000*eps)/binwidth1) + 1;
p1(p1>size(h_X2F_X2_X1,1)) = size(h_X2F_X2_X1,1);
p2(p2>size(h_X2F_X2_X1,1)) = size(h_X2F_X2_X1,1);
p3(p3>size(h_X2F_X2_X1,1)) = size(h_X2F_X2_X1,1);
h_X2F_X2_X1 = accumarray([p1 p2 p3],1,[nbin2-1 nbin2-1 nbin1-1]);   % 4-D histogram
h_X2F_X2_X1 = h_X2F_X2_X1 / sum(sum(sum(h_X2F_X2_X1)));      % normalization

% Calculating marginal histograms
h_X2_X1 = squeeze(sum(h_X2F_X2_X1,1));   % P(X2, X1)
h_X2F_X2 = squeeze(sum(h_X2F_X2_X1,3));   % P(X2F, X2)
h_X2 = squeeze(sum(h_X2_X1,2));          % P(X2)

% Calculating transfer entropy
% PaPb = h_X2F_X2_X1 .* repmat(h_X2',[nbin2-1 1 nbin1-1]);
% PcPd = permute(repmat(h_X2_X1,[1 1 nbin2-1]),[3 1 2]) .* repmat(h_X2F_X2,[1 1 nbin1-1]);
% pte = h_X2F_X2_X1 .* log2(PaPb./PcPd);
% pte = pte(:);
% TE = nansum(pte);

% Alternative calculation of transfer entropy (for bias correction purposes)
ph = h_X2F_X2 .* log2(h_X2F_X2);       % H(X2F, X2)
H_X2F_X2 = -nansum(ph(:));
ph = h_X2 .* log2(h_X2);       % H(X2F, X2)
H_X2 = -nansum(ph(:));
H_X2FcX2 = H_X2F_X2 - H_X2;     % H(X2F|X2)
ph = h_X2F_X2_X1 .* log2(h_X2F_X2_X1);    % H(X2F, X2, X1)
H_X2F_X2_X1 = -nansum(ph(:));
ph = h_X2_X1 .* log2(h_X2_X1);    % H(X2,X1)
H_X2_X1 = -nansum(ph(:));
H_X2FcX2_X1 = H_X2F_X2_X1 - H_X2_X1;    % H(X2F|X1,X2)
TE = H_X2FcX2 - H_X2FcX2_X1;

% Bias correction
Nt = n;   % total number of trials
Rs_bar = zeros(1,size(h_X2F_X2,2));
for k2 = 1:size(h_X2F_X2,2)
    Rs_bar(k2) = bayescount(Nt,h_X2F_X2(:,k2));
end
Bias_HRS = ((-1) / (2 * Nt * log(2))) * sum(Rs_bar-1);
H_X2FcX2_corr = H_X2FcX2 - Bias_HRS;

Nt = n;   % total number of trials
Rs_bar = zeros(size(h_X2F_X2_X1,2),size(h_X2F_X2_X1,3));
for k2 = 1:size(h_X2F_X2_X1,2)
    for k3 = 1:size(h_X2F_X2_X1,3)
        Rs_bar(k2,k3) = bayescount(Nt,h_X2F_X2_X1(:,k2,k3));
    end
end
Bias_HRS = ((-1) / (2 * Nt * log(2))) * sum(sum(Rs_bar-1));
H_X2FcX2_X1_corr = H_X2FcX2_X1 - Bias_HRS;

TE_corr = H_X2FcX2_corr - H_X2FcX2_X1_corr;
TE_bias = TE - TE_corr;



% -------------------------------------------------------------------------
function TE = ltren2(W1,W2,minW1,maxW1,minW2,maxW2,sr,tau)
%TREN   Transfer entropy.
%   TE = TREN(W1,W2,SR) calculates transfer entropy for time series W1 and
%   W2 sampled on SR (TE_W1->W2).
%
%   TE = TREN(W1,W2,SR,TAU) uses TAU ms as time lag between current and
%   "future" values.
%
%   [TE TE_CORR TE_BIAS] = TREN(W1,W2,SR,TAU) returns unbiased estimate of
%   transfer entropy applying Treves-Panzeri bias correction algorithm. The
%   amount of the bias of the original estimate is also returned.
%
%   [TE TE_CORR TE_BIAS H] = TREN(W1,W2,SR,TAU) returns H(X2F|X2)
%   conditonal entropy for transfer entropy normalization.
%
%   References:
%   (1) Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007)
%   Correcting for the sampling bias problem in spike train information
%   measures. J Neurophysiol 98:1064-1072.
%   (2) Gour�vitch B, Eggermont JJ (2007) Evaluating information transfer
%   between auditory cortical neurons. J Neurophysiol 97:2533-2543.
%   (3) Imas OA, Ropella KM, Ward BD, Wood JD, Hudetz AG (2005) Volatile
%   anesthetics disrupt frontal-posterior recurrent information transfer at
%   gamma frequencies in rat. Neurosci Lett 387:145-50.

% Input argument check
error(nargchk(7,8,nargin))
if nargin < 8
    tau = 100;   % time lag in ms
end

% Extracting the random variables
tau2 = tau / 1000 * sr;     % time lag in data points
X1 = W1(:,1:end-tau2);
X2 = W2(:,1:end-tau2);
X2F = W2(:,tau2+1:end);
n = numel(X1);

% Calculating joint histograms
bno = 30;   % bin number for histogram estimates
minX1 = minW1;    % discretization of X1
maxX1 = maxW1;
binwidth1 = (maxX1 - minX1) ./ bno;
xx1 = minX1 + binwidth1 * (0:bno);   % bin edges
xx1(length(xx1)) = maxX1;
xx1(1) = -inf;
nbin1 = length(xx1);

minX2 = minW2;    % discretization of X2
maxX2 = maxW2;
binwidth2 = (maxX2 - minX2) ./ bno;
xx2 = minX2 + binwidth2 * (0:bno);   % bin edges
xx2(length(xx2)) = maxX2;
xx2(1) = -inf;
nbin2 = length(xx2);

h_X2F_X2_X1 = zeros(nbin2-1,nbin2-1,nbin1-1);   % P(X2F, X2, X1)
t1 = X2F(:) - minX2;
t2 = X2(:) - minX2;
t3 = X1(:) - minX1;
p1 = fix((t1-1000000*eps)/binwidth2) + 1;
p2 = fix((t2-1000000*eps)/binwidth2) + 1;
p3 = fix((t3-1000000*eps)/binwidth1) + 1;
p1(p1>size(h_X2F_X2_X1,1)) = size(h_X2F_X2_X1,1);
p2(p2>size(h_X2F_X2_X1,1)) = size(h_X2F_X2_X1,1);
p3(p3>size(h_X2F_X2_X1,1)) = size(h_X2F_X2_X1,1);
h_X2F_X2_X1 = accumarray([p1 p2 p3],1,[nbin2-1 nbin2-1 nbin1-1]);   % 4-D histogram
h_X2F_X2_X1 = h_X2F_X2_X1 / sum(sum(sum(h_X2F_X2_X1)));      % normalization

% Calculating marginal histograms
h_X2_X1 = squeeze(sum(h_X2F_X2_X1,1));   % P(X2, X1)
h_X2F_X2 = squeeze(sum(h_X2F_X2_X1,3));   % P(X2F, X2)
h_X2 = squeeze(sum(h_X2_X1,2));          % P(X2)

% Calculating transfer entropy
PaPb = h_X2F_X2_X1 .* repmat(h_X2',[nbin2-1 1 nbin1-1]);
PcPd = permute(repmat(h_X2_X1,[1 1 nbin2-1]),[3 1 2]) .* repmat(h_X2F_X2,[1 1 nbin1-1]);
pte = h_X2F_X2_X1 .* log2(PaPb./PcPd);
pte = pte(:);
TE = nansum(pte);

% H(X2F|X2)
% ph = h_X2F_X2 .* log2(h_X2F_X2);       % H(X2F, X2)
% H_X2F_X2 = -nansum(ph(:));
% ph = h_X2 .* log2(h_X2);       % H(X2F, X2)
% H_X2 = -nansum(ph(:));
% H_X2FcX2 = H_X2F_X2 - H_X2;     % H(X2F|X2)