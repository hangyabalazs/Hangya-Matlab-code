function trtrenrun_noTP
%TRTRENRUN_NOTP   Calculates time-resolved transfer entropy.
%   TRTRENRUN_NOTP calculates wavelet power transfer entropy for the
%   longest theta and non-theta segments in one-second windows. (The
%   segments are attached together for wavelet calculation and extrema for
%   discretisation are calculated from the whole wavelet.) See TREN for
%   detailed calculations and references.
%
%   TRENRUN3_NOTP does not apply TP bias correction method.
%
%   See also TREN and TRENRUN.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
global DATADIR2
inpdir1 = [DATAPATH 'Burst\Cluster\Theta\'];   % input directory
inpdir2 = [DATAPATH 'Burst\Cluster\Noth\'];   % input directory
inpdir_tau = [DATAPATH 'TE\second_try_tauNTEmax_noTP\'];   % input directory
resdir = [DATAPATH 'TE\timeresolved_tauNTEmax_noTP\'];
% datadir = [DATADIR 'mi_zshift_data\discriminated_hc\'];
datadir = DATADIR2;
mm = pwd;
cd(resdir)

% Import
[files1 files1_short] = b_filelist(inpdir1);
[files2 files2_short] = b_filelist(inpdir2);
files_short = intersect(files1_short,files2_short);
% files_short = files_short(4)

sf = length(files_short);
[datalist, dlist_short] = b_filelist(datadir);

% Progress indicator
wb = waitbar(0,'Running TRTRENRUN...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Initialization
mDF_theta = zeros(1,sf);
mDF_noth = zeros(1,sf);

% Load
for o = 1:sf
    ftheta = files1(find(strcmp(files_short{o},files1_short))).name;
    fnoth = files2(find(strcmp(files_short{o},files2_short))).name;
    filenam = ftheta(1:6);
    if ~isequal(filenam,fnoth(1:6))
        error('Technical error 47')
    end
    cmps = strread(ftheta,'%s','delimiter','_MH');   % M, H: delimiters for 3ch data
    seglen_theta = str2num(cmps{5}) - str2num(cmps{4});
    i_first_theta = str2num(cmps{4});
    i_second_theta = str2num(cmps{5});
    cmps = strread(ftheta,'%s','delimiter','_');
    ftheta = cmps{1};
    for kt = 2:length(cmps)-1
        ftheta = [ftheta '_' cmps{kt}];
    end
    cmps = strread(fnoth,'%s','delimiter','_MH');   % M, H: delimiters for 3ch data
    seglen_noth = str2num(cmps{5}) - str2num(cmps{4});
    i_first_noth = str2num(cmps{4});
    i_second_noth = str2num(cmps{5});
    cmps = strread(fnoth,'%s','delimiter','_');
    fnoth = cmps{1};
    for kn = 2:length(cmps)-1
        fnoth = [fnoth '_' cmps{kn}];
    end
    
    inx = find(strcmp(filenam,dlist_short));
    fn = datalist(inx).name;
    ffn = [datadir fn];
    load(ffn);      % load raw data
    
    ff = [inpdir_tau filenam '_mTE.mat'];    % load optimal tau
    load(ff)
    
% Create new data (attach theta and non-theta segment)
    eeg_theta = eeg(i_first_theta:i_second_theta);
    eeg_noth = eeg(i_first_noth:i_second_noth);
    eeg_new = [eeg_noth eeg_theta];
    vdisc_theta = vdisc(find(vdisc>=i_first_theta&vdisc<=i_second_theta));
    vdisc_noth = vdisc(find(vdisc>=i_first_noth&vdisc<=i_second_noth));
    vdisc_new = [vdisc_noth-i_first_noth vdisc_theta-i_first_theta+i_second_noth-i_first_noth];
    seglen_new = seglen_theta + seglen_noth;
    
% Create random unit (shuffle ISIs)
    isi_theta = diff(vdisc_theta);
    isi_noth = diff(vdisc_noth);
    
    lit = length(isi_theta);
    rp = randperm(lit);
    while any(rp==(1:lit))
        rp = randperm(lit);
    end
    psi1 = [];
    for it = 1:lit
        psi1 = [psi1 isi_theta(rp(it))];
    end
    psvd_theta = [vdisc_theta(1) vdisc_theta(1)+cumsum(psi1)];
    
    lin = length(isi_noth);
    rp = randperm(lin);
    while any(rp==(1:lin))
        rp = randperm(lin);
    end
    psi2 = [];
    for in = 1:lin
        psi2 = [psi2 isi_noth(rp(in))];
    end
    psvd_noth = [vdisc_noth(1) vdisc_noth(1)+cumsum(psi2)];
    psvd_new = [psvd_noth-i_first_noth psvd_theta-i_first_theta+i_second_theta-i_first_theta];
    
% Downsample & wavelet transformation of eeg
    wavea_abs = eeg_wavelet(eeg_new(1:10:end)); % eeg (downsample on 1000 Hz)

% Sinc convolution & wavelet transformation of unit
    waveb_abs = unit_wavelet(vdisc_new,seglen_new);   % unit
    [wavec_abs,wavec_phase,f] = unit_wavelet(psvd_new,seglen_new);    % random unit
    
    if size(wavea_abs,2) > size(waveb_abs,2)    % adjust size
        wavea_abs = wavea_abs(:,1:end-1);
    elseif size(waveb_abs,2) > size(wavea_abs,2)
        waveb_abs = waveb_abs(:,1:end-1);
        wavec_abs = wavec_abs(:,1:end-1);
    end

% Transfer entropy
    fnd = find(f>6);    % frequency band bounderies
    pwind1 = fnd(end);
    fnd = find(f<2.5);
    pwind2 = fnd(1);
    
    abs1 = wavea_abs(pwind1:pwind2,:);  % eeg
    min_abs1 = min(abs1(:));
    max_abs1 = max(abs1(:));
    abs1_theta = abs1(:,round((i_second_noth-i_first_noth)/10)+1:end);
    abs1_noth = abs1(:,1:round((i_second_noth-i_first_noth)/10));
    abs2 = waveb_abs(pwind1:pwind2,:);  % unit
    min_abs2 = min(abs2(:));
    max_abs2 = max(abs2(:));
    abs2_theta = abs2(:,round((i_second_noth-i_first_noth)/10)+1:end);
    abs2_noth = abs2(:,1:round((i_second_noth-i_first_noth)/10));
    abs3 = wavec_abs(pwind1:pwind2,:);  % random unit
    min_abs3 = min(abs3(:));
    max_abs3 = max(abs3(:));
    abs3_theta = abs3(:,round((i_second_noth-i_first_noth)/10)+1:end);
    abs3_noth = abs3(:,1:round((i_second_noth-i_first_noth)/10));
    
    next = 1;
    [k1 k2] = size(abs1_theta);
    winlen = 1000;   % window size
    maxi = floor(k2/winlen);
    ovlp = 1;
    tau = mtau_theta;
    mmax = maxi * ovlp - ovlp + 1;
    TE_theta_eu = zeros(1,mmax);
    TE_theta_eu_shuffled = zeros(1,mmax);
    NTE_theta_eu = zeros(1,mmax);
    TE_theta_ue = zeros(1,mmax);
    TE_theta_ue_shuffled = zeros(1,mmax);
    NTE_theta_ue = zeros(1,mmax);
    for i = 1:mmax        % THETA LOOP
        inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
        inx1 = round(inx1);
        inx2 = inx1 + winlen - 1;
        [TE H_X2FcX2] = ltren(abs1_theta(:,inx1:inx2),abs2_theta(:,inx1:inx2),...
            min_abs1,max_abs1,min_abs2,max_abs2,1000,tau);  % theta, EEG->unit
        TE_theta_eu(i) = TE;
        TE = ltren(abs1_theta(:,inx1:inx2),abs3_theta(:,inx1:inx2),...
            min_abs1,max_abs1,min_abs3,max_abs3,1000,tau);
        TE_theta_eu_shuffled(i) = TE;
        NTE_theta_eu(i) = (TE_theta_eu(i) - TE_theta_eu_shuffled(i)) / H_X2FcX2;
        NTE_theta_eu(i) = max(NTE_theta_eu(i),0);
        
        [TE H_X2FcX2] = ltren(abs2_theta(:,inx1:inx2),abs1_theta(:,inx1:inx2),...
            min_abs2,max_abs2,min_abs1,max_abs1,1000,tau);  % theta, unit->EEG
        TE_theta_ue(i) = TE;
        TE = ltren(abs3_theta(:,inx1:inx2),abs1_theta(:,inx1:inx2),...
            min_abs3,max_abs3,min_abs1,max_abs1,1000,tau);
        TE_theta_ue_shuffled(i) = TE;
        NTE_theta_ue(i) = (TE_theta_ue(i) - TE_theta_ue_shuffled(i)) / H_X2FcX2;
        NTE_theta_ue(i) = max(NTE_theta_ue(i),0);
        
        next = next + 1;
    end
    DF_theta_eu = (NTE_theta_eu - NTE_theta_ue) ./ (NTE_theta_eu + NTE_theta_ue);    % direction of flow (pos: EEG->unit)
    DIF_theta_eu = (NTE_theta_eu - NTE_theta_ue);
    
    [k1 k2] = size(abs1_noth);
    maxi = floor(k2/winlen);
    tau = mtau_noth;
    mmax = maxi * ovlp - ovlp + 1;
    TE_noth_eu = zeros(1,mmax);
    TE_noth_eu_shuffled = zeros(1,mmax);
    NTE_noth_eu = zeros(1,mmax);
    TE_noth_ue = zeros(1,mmax);
    TE_noth_ue_shuffled = zeros(1,mmax);
    NTE_noth_ue = zeros(1,mmax);
    for i = 1:mmax        % NON-THETA LOOP
        inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
        inx1 = round(inx1);
        inx2 = inx1 + winlen - 1;
        [TE H_X2FcX2] = ltren(abs1_noth(:,inx1:inx2),abs2_noth(:,inx1:inx2),...
            min_abs1,max_abs1,min_abs2,max_abs2,1000,tau);  % non-theta, EEG->unit
        TE_noth_eu(i) = TE;
        TE = ltren(abs1_noth(:,inx1:inx2),abs3_noth(:,inx1:inx2),...
            min_abs1,max_abs1,min_abs3,max_abs3,1000,tau);
        TE_noth_eu_shuffled(i) = TE;
        NTE_noth_eu(i) = (TE_noth_eu(i) - TE_noth_eu_shuffled(i)) / H_X2FcX2;
        NTE_noth_eu(i) = max(NTE_noth_eu(i),0);
        
        [TE H_X2FcX2] = ltren(abs2_noth(:,inx1:inx2),abs1_noth(:,inx1:inx2),...
            min_abs2,max_abs2,min_abs1,max_abs1,1000,tau);  % non-theta, unit->EEG
        TE_noth_ue(i) = TE;
        TE = ltren(abs3_noth(:,inx1:inx2),abs1_noth(:,inx1:inx2),...
            min_abs3,max_abs3,min_abs1,max_abs1,1000,tau);
        TE_noth_ue_shuffled(i) = TE;
        NTE_noth_ue(i) = (TE_noth_ue(i) - TE_noth_ue_shuffled(i)) / H_X2FcX2;
        NTE_noth_ue(i) = max(NTE_noth_ue(i),0);
        
        next = next + 1;
    end
    DF_noth_eu = (NTE_noth_eu - NTE_noth_ue) ./ (NTE_noth_eu + NTE_noth_ue);
    DIF_noth_eu = (NTE_noth_eu - NTE_noth_ue);
    
% Plot and save
    H1 = figure;
    plot(DF_theta_eu,'r')
    title('DF theta EEG->unit')
    ffn = [filenam '_DFtheta.fig'];
    saveas(H1,ffn)
    H2 = figure;
    plot(DF_noth_eu)
    title('DF non-theta EEG->unit')
    ffn = [filenam '_DFnoth.fig'];
    saveas(H2,ffn)
    H1b = figure;
    plot(DIF_theta_eu,'r')
    title('DIF theta EEG->unit')
    ffn = [filenam '_DIFtheta.fig'];
    saveas(H1b,ffn)
    H2b = figure;
    plot(DIF_noth_eu)
    title('DIF non-theta EEG->unit')
    ffn = [filenam '_DIFnoth.fig'];
    saveas(H2b,ffn)
    mfn = [filenam '_TE.mat'];
    save(mfn,'DF_theta_eu','DF_noth_eu','TE_theta_eu','TE_theta_ue','TE_noth_eu','TE_noth_ue',...
        'NTE_theta_eu','NTE_theta_ue','NTE_noth_eu','NTE_noth_ue');
    close all
    clear wavea_abs waveb_abs wavec_abs r s
    waitbar(o/sf)
end
% H3 = figure;
% hold on
% plot(mTau_noth,mDF_noth,'b.')
% plot(mTau_theta,mDF_theta,'r.')
% xlabel('Tau max')
% ylabel('DF EEG->unit max')
% ffn = 'mDF.fig';
% saveas(H3,ffn)
% save('mDF.mat','mTau_theta','mTau_noth','mDF_theta','mDF_noth')

close(wb)
clear global FTHETA FNOTH DATINX1_THETA DATINX2_THETA DATINX1_NOTH DATINX2_NOTH
cd(mm)



% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function [pow,phase,f] = eeg_wavelet(dat)

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
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);



% -------------------------------------------------------------------------
function [pow,phase,f] = unit_wavelet(vdisc,lenu)

% Sinc convolution
fs = 10000;     % unit
dto = 1 / fs;
ts = zeros(1,lenu);
ts(vdisc) = 1;
du = diff(vdisc);
fdu = 1 ./ du;
fdu = [fdu 0.0001];
fcut = 100; 
fsnew = 1000;
dtnew = 1 / 1000;
fsold = 10000;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end

% Prepare for wavelet transformation
variance = std(zint) ^ 2;
zint = (zint - mean(zint)) / sqrt(variance) ;
n = length(zint);
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
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(zint,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);



% -------------------------------------------------------------------------
% ENTROPY
% -------------------------------------------------------------------------
function [TE H_X2FcX2] = ltren(W1,W2,minW1,maxW1,minW2,maxW2,sr,tau)
%TREN   Transfer entropy.
%   TE = TREN(W1,W2,SR) calculates transfer entropy for time series W1 and
%   W2 sampled on SR (TE_W1->W2).
%
%   TE = TREN(W1,W2,SR,TAU) uses TAU ms as time lag between current and
%   "future" values.
%
%   [TE H] = TREN(W1,W2,SR,TAU) returns H(X2F|X2)
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
for k = 1:n
    ind1 = min(p1(k),size(h_X2F_X2_X1,1));
    ind2 = min(p2(k),size(h_X2F_X2_X1,2));
    ind3 = min(p3(k),size(h_X2F_X2_X1,3));
    h_X2F_X2_X1(ind1,ind2,ind3) = h_X2F_X2_X1(ind1,ind2,ind3) + 1;
end
h_X2F_X2_X1 = h_X2F_X2_X1 / sum(sum(sum(h_X2F_X2_X1)));      % normalization

% Calculating marginal histograms
h_X2_X1 = squeeze(sum(h_X2F_X2_X1,1));   % P(X2, X1)
h_X2F_X2 = squeeze(sum(h_X2F_X2_X1,3));   % P(X2F, X2)
h_X2 = squeeze(sum(h_X2_X1,2));          % P(X2)

% Calculating transfer entropy
TE = 0;
for k1 = 1:nbin2-1
    for k2 = 1:nbin2-1
        for k3 = 1:nbin1-1
            if h_X2F_X2_X1(k1,k2,k3) ~= 0
                Pa = h_X2F_X2_X1(k1,k2,k3);
                Pb = h_X2(k2);
                Pc = h_X2_X1(k2,k3);
                Pd = h_X2F_X2(k1,k2);
                TE = TE + Pa * log2((Pa*Pb)/(Pc*Pd));
            end
        end
    end
end

% Alternative calculation of transfer entropy (for bias correction purposes)
H_X2F_X2 = 0;       % H(X2F, X2)
for k1 = 1:nbin2-1
    for k2 = 1:nbin2-1
        if h_X2F_X2(k1,k2) ~= 0
            Pa = h_X2F_X2(k1,k2);
            H_X2F_X2 = H_X2F_X2 - Pa * log2(Pa);
        end
    end
end
H_X2 = 0;           % H(X2)
for k2 = 1:nbin2-1
    if h_X2(k2) ~= 0
        Pa = h_X2(k2);
        H_X2 = H_X2 - Pa * log2(Pa);
    end
end
H_X2FcX2 = H_X2F_X2 - H_X2;     % H(X2F|X2)

H_X2F_X2_X1 = 0;    % H(X2F, X2, X1)
for k1 = 1:nbin2-1
    for k2 = 1:nbin2-1
        for k3 = 1:nbin1-1
            if h_X2F_X2_X1(k1,k2,k3) ~= 0
                Pa = h_X2F_X2_X1(k1,k2,k3);
                H_X2F_X2_X1 = H_X2F_X2_X1 - Pa * log2(Pa);
            end
        end
    end
end
H_X2_X1 = 0;    % H(X2,X1)
for k2 = 1:nbin2-1
    for k3 = 1:nbin1-1
        if h_X2_X1(k2,k3) ~= 0
            Pa = h_X2_X1(k2,k3);
            H_X2_X1 = H_X2_X1 - Pa * log2(Pa);
        end
    end
end
H_X2FcX2_X1 = H_X2F_X2_X1 - H_X2_X1;    % H(X2F|X1,X2)
TE2 = H_X2FcX2 - H_X2FcX2_X1;