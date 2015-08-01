function mspoptimaltau
%MSPOPTIMALTAU   Optimal time lag for TE calculation.
%   MSPOPTIMALTAU calculates transfers entropy for the total length of
%   longest theta and nontheta segments. It saves optimal tau values for
%   time-resolved TE calculations.
%
%   See also MSPTRTE.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
ratno = '4';
inpdir = [DATADIR 'MSHCsp\Viktor' ratno '\'];
thetadir = [DATAPATH 'MSHCsp\Wavelet\theta_segments\'];
nonthetadir = [DATAPATH 'MSHCsp\Wavelet\nontheta_segments\'];
resdir = [DATAPATH 'MSHCsp\TE\Viktor' ratno '\'];
mm = pwd;

% Filelist
switch  ratno
    case '4'
        rid = '1234';
        flist = {'201003171' '201003172' '201003173' '201003174'};
    case '5'
        rid = '1237'; %Viktor5's files
        flist = {'201003291' '201003292' '201003293' '201003294'};
end
sf = length(flist);

% Main
sr = 20000;    % sampling rate
dsr = 1000;    % resample at 1000 Hz
const = sr / dsr;
Inx = 33:65;
Inx = 65;
for o = 3:sf
    
    % Load EEG
    fname = flist{o};
    ff = [thetadir 'THETA_SEGMENTS_' rid '_' fname '.mat'];     % load theta segments
    load(ff)
    ThetaSegments = uniteseg(ThetaSegments,sr);     % drop gaps < 0.5 s
    ThetaSegments = short_killer(ThetaSegments);    % drop segments < 3 s
    ThetaSegments =  round(ThetaSegments/const);    % change sampling rate to 'dsr'
    lent = ThetaSegments(2,:) - ThetaSegments(1,:);     % theta segment length
    minx = find(lent==max(lent));   % find the longest theta segment
    th1 = ThetaSegments(1,minx);
    th2 = ThetaSegments(2,minx);

    ff = [nonthetadir 'NONTHETA_SEGMENTS_' rid '_' fname '.mat'];     % load non-theta segments
    load(ff)
    NonThetaSegments = short_killer(NonThetaSegments);    % drop segments < 3 s
    NonThetaSegments =  round(NonThetaSegments/const);    % change sampling rate to 'dsr'
    lenn = NonThetaSegments(2,:) - NonThetaSegments(1,:);     % length of non-theta segments
    minx = find(lenn==max(lenn));   % find the longest non-theta segment
    nth1 = NonThetaSegments(1,minx);
    nth2 = NonThetaSegments(2,minx);
    
    ff = [inpdir flist{o} '.par'];  % load channel settings
    pars = LoadPar(ff);
    pchs = cellfun(@(x) x',pars.ElecGp,'UniformOutput',0);
    chs = cell2mat(pchs);
    
    fnam = [inpdir flist{o} '.dat'];    % load EEG
    numch = pars.nChannels;
    buffersize = 2^14;
    fileinfo = dir(fnam);  % get file size and calculate the number of samples per channel
    flen = ceil(fileinfo(1).bytes/2/numch);
    numelmax = flen;
%     numelmin = flen - 1000000;
    numelmin = 0;
    datafile = fopen(fnam,'r');
    eeg = zeros(numelmax-numelmin+1,1);

    for ch = Inx % get the eeg channel

        numel = 0;
        fseek(datafile,numelmin*numch*2,'bof');
        while numelmin + numel < numelmax
            [data,count] = fread(datafile,[numch,buffersize],'int16');
            numelm = count / numch;
            chselect = chs(ch) + 1;
            eeg(numel+1:numel+numelm) = (data(chselect,:))';

            numel = numel + numelm;
        end

%         ff = [inpdir 'EEG\' rid '_' flist{o} '_eeg.mat'];
%         load(ff)
%         ch = 60;

        % Resample EEG
        eeg = eeg(1:const:size(eeg,1));     % downsample at 'dsr'

        % Load unit
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

            for nc = 2:ncells
                vdisc = res(clu==nc)';
                vdisc = round(vdisc/const);  % resample at 'dsr'
                vdisc_theta = vdisc(vdisc>th1&vdisc<th2);
                vdisc_noth = vdisc(vdisc>nth1&vdisc<nth2);

                % TRANSFER ENTROPY
                % Create new data (attach theta and non-theta segment)
                eeg_theta = eeg(th1:th2)';
                eeg_noth = eeg(nth1:nth2)';
                eeg_new = [eeg_noth eeg_theta];
                vdisc_new = [vdisc_noth-nth1 vdisc_theta-th1+nth2-nth1];
                seglen_theta = th2 - th1;
                seglen_noth = nth2 - nth1;
                seglen_new = seglen_theta + seglen_noth;

                % Create random unit (shuffle ISIs)
                [psvd_theta psvd_noth] = isi_shuffle(vdisc_theta,vdisc_noth);
                psvd_new = [psvd_noth-nth1 psvd_theta-th1+nth2-nth1];

                % Wavelet transformation of the EEG
                abs1 = eeg_wavelet(eeg_new);

                % Sinc convolution & wavelet transformation of unit
                abs2 = unit_wavelet(vdisc_new,seglen_new,dsr,dsr);   % unit
                abs3 = unit_wavelet(psvd_new,seglen_new,dsr,dsr);    % random unit

                if ~isequal(size(abs1,2),size(abs2,2))    % adjust size
                    mininx = min([size(abs1,2) size(abs2,2)]);
                    abs1 = abs1(:,1:mininx);
                    abs2 = abs2(:,1:mininx);
                    abs3 = abs3(:,1:mininx);
                end

                % Transfer entropy
                min_abs1 = min(abs1(:));    % calculate extrema
                max_abs1 = max(abs1(:));
                abs1_theta = abs1(:,nth2-nth1+1:end);
                abs1_noth = abs1(:,1:nth2-nth1);
                min_abs2 = min(abs2(:));
                max_abs2 = max(abs2(:));
                abs2_theta = abs2(:,nth2-nth1+1:end);
                abs2_noth = abs2(:,1:nth2-nth1);
                min_abs3 = min(abs3(:));
                max_abs3 = max(abs3(:));
                abs3_theta = abs3(:,nth2-nth1+1:end);
                abs3_noth = abs3(:,1:nth2-nth1);

                T = 0:10:200;
                T(1) = 1;   % preallocate 'current' type variables
                cTE = struct('theta',struct('eu',[],'ue',[]),'noth',struct('eu',[],'ue',[]));
                cNTE = struct('theta',struct('eu',[],'ue',[]),'noth',struct('eu',[],'ue',[]));
                cBias = struct('theta',struct('eu',[],'eu_shuffled',[],'ue',[],'ue_shuffled',[]),...
                    'noth',struct('eu',[],'eu_shuffled',[],'ue',[],'ue_shuffled',[]));
                cDF_theta_eu = zeros(1,length(T));
                cDF_noth_eu = zeros(1,length(T));
                cNTE_theta_eu = zeros(1,length(T));
                cNTE_theta_ue = zeros(1,length(T));
                cNTE_noth_eu = zeros(1,length(T));
                cNTE_noth_ue = zeros(1,length(T));
                next = 1;
                for tau = T
%                     disp(tau)
                    [TE TE_corr TE_bias H_X2FcX2] = ltren(abs1_theta,abs2_theta,...
                        min_abs1,max_abs1,min_abs2,max_abs2,1000,tau);  % theta, EEG->unit
                    TE_theta_eu = TE;
                    TE_theta_eu_bias = TE_bias;
                    [TE TE_corr TE_bias] = ltren(abs1_theta,abs3_theta,...
                        min_abs1,max_abs1,min_abs3,max_abs3,1000,tau);
                    TE_theta_eu_shuffled = TE;
                    TE_theta_eu_shuffled_bias = TE_bias;
                    NTE_theta_eu = (TE_theta_eu - TE_theta_eu_shuffled) / H_X2FcX2;
                    NTE_theta_eu = max(NTE_theta_eu,0);

                    [TE TE_corr TE_bias H_X2FcX2] = ltren(abs2_theta,abs1_theta,...
                        min_abs2,max_abs2,min_abs1,max_abs1,1000,tau);  % theta, unit->EEG
                    TE_theta_ue = TE;
                    TE_theta_ue_bias = TE_bias;
                    [TE TE_corr TE_bias] = ltren(abs3_theta,abs1_theta,...
                        min_abs3,max_abs3,min_abs1,max_abs1,1000,tau);
                    TE_theta_ue_shuffled = TE;
                    TE_theta_ue_shuffled_bias = TE_bias;
                    NTE_theta_ue = (TE_theta_ue - TE_theta_ue_shuffled) / H_X2FcX2;
                    NTE_theta_ue = max(NTE_theta_ue,0);

                    [TE TE_corr TE_bias H_X2FcX2] = ltren(abs1_noth,abs2_noth,...
                        min_abs1,max_abs1,min_abs2,max_abs2,1000,tau);  % non-theta, EEG->unit
                    TE_noth_eu = TE;
                    TE_noth_eu_bias = TE_bias;
                    [TE TE_corr TE_bias] = ltren(abs1_noth,abs3_noth,...
                        min_abs1,max_abs1,min_abs3,max_abs3,1000,tau);
                    TE_noth_eu_shuffled = TE;
                    TE_noth_eu_shuffled_bias = TE_bias;
                    NTE_noth_eu = (TE_noth_eu - TE_noth_eu_shuffled) / H_X2FcX2;
                    NTE_noth_eu = max(NTE_noth_eu,0);

                    [TE TE_corr TE_bias H_X2FcX2] = ltren(abs2_noth,abs1_noth,...
                        min_abs2,max_abs2,min_abs1,max_abs1,1000,tau);  % non-theta, unit->EEG
                    TE_noth_ue = TE;
                    TE_noth_ue_bias = TE_bias;
                    [TE TE_corr TE_bias] = ltren(abs3_noth,abs1_noth,...
                        min_abs3,max_abs3,min_abs1,max_abs1,1000,tau);
                    TE_noth_ue_shuffled = TE;
                    TE_noth_ue_shuffled_bias = TE_bias;
                    NTE_noth_ue = (TE_noth_ue - TE_noth_ue_shuffled) / H_X2FcX2;
                    NTE_noth_ue = max(NTE_noth_ue,0);

                    if abs(NTE_theta_ue) + abs(NTE_theta_eu) + abs(NTE_noth_ue) + abs(NTE_noth_eu) == Inf
                        keyboard
                    end

                    DF_theta_eu = (NTE_theta_eu - NTE_theta_ue) / (NTE_theta_eu + NTE_theta_ue);    % direction of flow (pos: EEG->unit)
                    DF_noth_eu = (NTE_noth_eu - NTE_noth_ue) / (NTE_noth_eu + NTE_noth_ue);

                    cTE.theta.eu = [cTE.theta.eu TE_theta_eu];
                    cTE.theta.ue = [cTE.theta.ue TE_theta_ue];
                    cTE.noth.eu = [cTE.noth.eu TE_noth_eu];
                    cTE.noth.ue = [cTE.noth.ue TE_noth_ue];
                    cNTE.theta.eu = [cNTE.theta.eu NTE_theta_eu];
                    cNTE.theta.ue = [cNTE.theta.ue NTE_theta_ue];
                    cNTE.noth.eu = [cNTE.noth.eu NTE_noth_eu];
                    cNTE.noth.ue = [cNTE.noth.ue NTE_noth_ue];
                    cBias.theta.eu = [cBias.theta.eu TE_theta_eu_bias];
                    cBias.theta.eu_shuffled = [cBias.theta.eu_shuffled TE_theta_eu_shuffled_bias];
                    cBias.theta.ue = [cBias.theta.ue TE_theta_ue_bias];
                    cBias.theta.ue_shuffled = [cBias.theta.ue_shuffled TE_theta_ue_shuffled_bias];
                    cBias.noth.eu = [cBias.noth.eu TE_noth_eu_bias];
                    cBias.noth.eu_shuffled = [cBias.noth.eu_shuffled TE_noth_eu_shuffled_bias];
                    cBias.noth.ue = [cBias.noth.ue TE_noth_ue_bias];
                    cBias.noth.ue_shuffled = [cBias.noth.ue_shuffled TE_noth_ue_shuffled_bias];
                    cDF_theta_eu(next) = DF_theta_eu;
                    cDF_noth_eu(next) = DF_noth_eu;
                    cNTE_theta_eu(next) = NTE_theta_eu;
                    cNTE_theta_ue(next) = NTE_theta_ue;
                    cNTE_noth_eu(next) = NTE_noth_eu;
                    cNTE_noth_ue(next) = NTE_noth_ue;
                    next = next + 1;
                end
                mab = max(abs(cNTE_theta_eu),abs(cNTE_theta_ue));
                thetainx = find(mab==max(mab));
                mab = max(abs(cNTE_noth_eu),abs(cNTE_noth_ue));
                nothinx = find(mab==max(mab));
                if length(thetainx) > 1
                    thetainx = thetainx(1);
                end
                if length(nothinx) > 1
                    nothinx = nothinx(1);
                end
                mdf_theta = cDF_theta_eu(thetainx);
                mDF_theta(o) = mdf_theta;
                mdf_noth = cDF_noth_eu(nothinx);
                mDF_noth(o) = mdf_noth;
                mtau_theta = T(thetainx);
                mTau_theta(o) = mtau_theta;
                mtau_noth = T(nothinx);
                mTau_noth(o) = mtau_noth;

                % Plot and save
                H1 = figure;
                plot(T,cDF_theta_eu,'r')
                title('DF theta EEG->unit')
                ffn = [resdir rid '_' fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc) '_DFtheta.fig' ];
                saveas(H1,ffn)
                H2 = figure;
                plot(T,cDF_noth_eu)
                title('DF non-theta EEG->unit')
                ffn = [resdir rid '_' fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc) '_DFnoth.fig' ];
                saveas(H2,ffn)
                mfn = [resdir rid '_' fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc) '_TE.mat' ];
                save(mfn,'cDF_theta_eu','cDF_noth_eu','cTE','cNTE','cBias')
                mfn = [resdir rid '_' fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc) '_mTE.mat' ];
                save(mfn,'mdf_theta','mdf_noth','mtau_theta','mtau_noth')
                close all
                clear wavea_abs waveb_abs wavec_abs

                % TE significance values
                snm = 100;
                sTE_theta_eu = zeros(1,snm);
                sTE_theta_ue = zeros(1,snm);
                sTE_noth_eu = zeros(1,snm);
                sTE_noth_ue = zeros(1,snm);
                for k = 1:snm
%                     disp(k)
                    [psvd_theta psvd_noth] = isi_shuffle(vdisc_theta,vdisc_noth); % create random unit (shuffle ISIs)
                    psvd_new = [psvd_noth-nth1 psvd_theta-th1+nth2-nth1];
                    abs4 = unit_wavelet(psvd_new,seglen_new,dsr,dsr); % sinc convolution & wavelet transformation of random unit
                    min_abs4 = min(abs4(:));
                    max_abs4 = max(abs4(:));
                    abs4_theta = abs4(:,nth2-nth1+1:end);
                    abs4_noth = abs4(:,1:nth2-nth1);

                    TE = ltren2(abs1_theta,abs4_theta,...
                        min_abs1,max_abs1,min_abs4,max_abs4,1000,mtau_theta);  % theta, EEG->unit
                    sTE_theta_eu(k) = TE;

                    TE = ltren2(abs4_theta,abs1_theta,...
                        min_abs4,max_abs4,min_abs1,max_abs1,1000,mtau_theta);  % theta, unit->EEG
                    sTE_theta_ue(k) = TE;

                    TE = ltren2(abs1_noth,abs4_noth,...
                        min_abs1,max_abs1,min_abs4,max_abs4,1000,mtau_noth);  % non-theta, EEG->unit
                    sTE_noth_eu(k) = TE;

                    TE = ltren2(abs4_noth,abs1_noth,...
                        min_abs4,max_abs4,min_abs1,max_abs1,1000,mtau_noth);  % non-theta, unit->EEG
                    sTE_noth_ue(k) = TE;
                end
                mfn = [resdir rid '_' fname '_' num2str(ch) '_' num2str(shankno) '_' num2str(nc) '_TEsign.mat' ];
                save(mfn,'sTE_theta_eu','sTE_theta_ue','sTE_noth_eu','sTE_noth_ue')
            end
        end
    end
end
fclose(datafile)
cd(mm)


% -------------------------------------------------------------------------
% SEGMENT SELECTION
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
% RANDOM UNIT
% -------------------------------------------------------------------------
function [psvd_theta psvd_noth] = isi_shuffle(vdisc_theta,vdisc_noth)

% Interspike intervals
isi_theta = diff(vdisc_theta);
isi_noth = diff(vdisc_noth);

% Randomize theta segment
lit = length(isi_theta);    % shuffle theta ISIs
rp = randperm(lit);
while any(rp==(1:lit))
    rp = randperm(lit);
end
psi1 = [];
for it = 1:lit
    psi1 = [psi1 isi_theta(rp(it))];
end
psvd_theta = [vdisc_theta(1) vdisc_theta(1)+cumsum(psi1)];

% Randomize non-theta segment
lin = length(isi_noth);     % shuffle non-theta ISIs
rp = randperm(lin);
while any(rp==(1:lin))
    rp = randperm(lin);
end
psi2 = [];
for in = 1:lin
    psi2 = [psi2 isi_noth(rp(in))];
end
psvd_noth = [vdisc_noth(1) vdisc_noth(1)+cumsum(psi2)];



% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function pow = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 1000;
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
function pow = unit_wavelet(vdisc,lenu,sr,dsr)

% Sinc convolution
fs = sr;     % unit
dto = 1 / fs;
fcut = 100; 
fsnew = dsr;
dtnew = 1 / fsnew;
fsratio = fsnew / fs;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);

tnew2 = (-lenu*fsratio:lenu*fsratio) * dtnew * fcut;
psinc = sinc(tnew2);
lt = length(tnew);
told2 = round(told/dtnew/fcut);
zint = zeros(1,lt);
for k = 1:lentold
    zint = zint + psinc(lt-told2(k)+2:2*lt-told2(k)+1);
end

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
	[daughter,fourier_factor,coi,dofmin] = wave_bases(mother,k,scale(a1),param);	
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
%   (2) Gourévitch B, Eggermont JJ (2007) Evaluating information transfer
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
%   (2) Gourévitch B, Eggermont JJ (2007) Evaluating information transfer
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