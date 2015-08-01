function noseg5
%NOSEG5   Analysis of theta and non-theta segments (NO project).
%   NOSEG5 calculates the wavelet power for non-theta centers in the last
%   1800 s for pre-injection files or in the first 3600 s in post-injection
%   files. Dirac functions centered to non-theta centers are sinc-convolved
%   and wavelet spectrum is calculated in the low frequency range to reveal
%   the rhytmic occurrence of non-theta segments. The results are plotted
%   against time.
%
%   See also NOSEG6.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
inpdir2 = [DATAPATH 'NO\Wavelet\Segments\'];
resdir = [DATAPATH 'NO\'];
tblfile = [DATAPATH 'NO\seg_data.xls'];

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSEG5...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Load segment data
[tbl0 tbl] = xlsread(tblfile);

% Main
sr = 5000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
newstep = 500;
dsr = sr / newstep;    % downsample on 10 Hz
names = cell(1,sf);
tps = zeros(1,sf);
for o = 3:sf
    fname = files(o).name;      % load data
    names{o} = fname;
    cmps = strread(fname(1:end-4),'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    
    load([inpdir fname])
    eeg = hEEG.values;
    len = length(eeg);
    
    fn2 = [fname(1:end-4) '_SEGMENTS'];     % load segment identification
    load([inpdir2 fn2])
    
    inx = find(strcmp({tbl{:,1}},fname));
    id1 = tbl{inx,2};
    id2 = tbl{inx,3};
    
    ntf = NonThetaSegments(1,:);    % wavelet for non-theta centers
    ntl = NonThetaSegments(2,:);
    ntc = (ntf + ntl) / 2;
    ntlen = (ntl - ntf) / sr;
    if strcmp(id2,'pre')
        ntc = ntc(ntc>len-1800*sr);
        lenu = 1800 * sr;
    elseif strcmp(id2,'post')
        ntc = ntc(ntc<3600*sr);
        lenu = 3600 * sr;
    end
    [powntc,f] = unitwavelet(round(ntc/newstep),round(min(lenu,len)/newstep),dsr);    % EEG wavelet
    H = figure;
    fx = find(f<0.005,1,'first');
%     imagesc(powntc./repmat(sum(powntc(fx:end,:)),size(powntc,1),1))
    imagesc(powntc)
%     set(gca,'CLim',[0 300])
    b_rescaleaxis('Y',f)
    title(titlestr)
    fnw = [resdir 'NtcWavelet\' fname(1:end-4) '_NTCWAVELET.jpg'];
    saveas(H,fnw)
    clear powntc
    close(H)
    waitbar(o/sf)
end
close(wb)

% Save
% xlswrite([resdir 'noallpow_pre1800post3600.xls'],names','sheet1','A1')
% xlswrite([resdir 'noallpow_pre1800post3600.xls'],tps','sheet1','B1')

% cd(mm)

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [pow,f] = unitwavelet(vdisc,lenu,sr)
%UNITWAVELET   Wavelet calculation.
%   [POW,PHASE,F] = UNITWAVELET(DATA,LENU,SR) performs wavelet calculation 
%   on discriminated unit data (VDISC) sampled at SR sampling rate for a
%   data length of LENU. Wavelet power (POW), phase (PHASE) and scale
%   vector (F) are returned. DATA is first convolved by a sinc kernel and
%   standardized. A minimal interesting frequency of 0.5 Hz is set.
%
%   See also EEGWAVELET.

% Sinc convolution
fs = sr;     % unit
dto = 1 / fs;
ts = zeros(1,lenu);
ts(vdisc) = 1;
du = diff(vdisc);
fdu = 1 ./ du;
fdu = [fdu 1/sr];
fcut = 100; 
fsnew = fs;
dtnew = 1 / fsnew;
fsold = sr;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc((tnew-told(i))/50);
end
% zint = gaussconv(vdisc,lenu);

% Prepare for wavelet transformation
variance = std(zint) ^ 2;
zint = (zint - mean(zint)) / sqrt(variance) ;
n = length(zint);
dt = 1 / fs;
pad = 1;
% spf = 150;
% omega0 = 6;
% c = 4 * pi / (omega0 + sqrt(2+omega0^2));
% scs = [1:-0.005:0.005];
% s = 1 ./ scs ./ c;
% fperiod = c .* s;
% f = 1 ./ fperiod;
dj = 0.08;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
s = s(70:end);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
% [wave,period,scale,coi] = wavelet(zint,dt,pad,s,mother,param);
[wave,period,scale,coi] = b_wavelet_new3(zint,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;

% -------------------------------------------------------------------------
function zint = gaussconv(vdisc,lenu)

ipunit = zeros(1,lenu);
ipunit(vdisc) = 1;
wbh = gausswin(200);
wipunit = conv(ipunit,wbh);
zint = wipunit(1:lenu);

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
function b_rescaleaxis(axis,scale)
%RESCALEAXIS   Changes ticklabeles to 'frequency'.
%   RESCALEAXIS(AXIS,SCALE) requires input parameter SCALE containing the correspondig
%   frequency value for each index. It changes the ticklabels of AXIS to the
%   frequency values.
%
%   See also RETIMEFIG and RESCALEFIG.

%Input arguments check
error(nargchk(2,2,nargin));

% Get label
if strcmp(axis,'X')
    a0 = get(gca,'XTick')';
elseif strcmp(axis,'Y')
    a0 = get(gca,'YTick')';
else
    error('First input argument must be ''X'' or ''Y''.')
end

% Interpolation
lo0 = floor(a0);
if ~isequal(a0,lo0);
    hi0 = ceil(a0);
    iv = hi0 - lo0;
    liv = a0 - lo0;
    fiv = find(iv);
    r = zeros(size(iv));
    r(fiv) = liv(fiv) ./ iv(fiv);
    milo = lo0 < 1 | lo0 > length(scale);
    mihi = hi0 < 1 | hi0 > length(scale);
    lo1 = zeros(size(lo0));
    hi1 = zeros(size(hi0));
    lo1(find(milo)) = NaN;
    hi1(find(mihi)) = NaN;
    lo1(find(~milo)) = scale(lo0(find(~milo)));
    hi1(find(~mihi)) = scale(hi0(find(~mihi)));
    a1 = lo1 + r .* (hi1 - lo1);
else
    mia = a0 < 1 | a0 > length(scale);
    a1 = zeros(size(a0));
    a1(find(mia)) = NaN;
    a1(find(~mia)) = scale(a0(find(~mia)));
end
a1 = round(a1*10000) / 10000;

% Set label
na1 = num2str(a1);
[m n] = size(na1);
rna = reshape(na1',1,m*n);
ni = findstr(rna,'NaN');
rna(ni) = ' ';
rna(ni+1) = ' ';
rna(ni+2) = ' ';
na1 = reshape(rna,n,m)';
if strcmp(axis,'X')
    set(gca,'XTickLabel',na1)
elseif strcmp(axis,'Y')
    set(gca,'YTickLabel',na1)
end