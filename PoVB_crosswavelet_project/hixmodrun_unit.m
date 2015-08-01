function hixmodrun_unit
%HIXMODRUN_UNIT   Runs cross-modulation on a sequence of files.
%   HIXMODRUN_UNIT calculates and saves crossmodulation plots (see Tort et al.,
%   2008, PNAS).
%
%   See also WAVELET_NEW3.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = [DATAPATH 'Hajni\EEGMPO\mat\'];
resdir = [DATAPATH 'Hajni\EEGMPO_new\crossmod_unit\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running HIXMODRUN UNIT...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
main(inpdir,resdir,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir,resdir,files_short,sf,wb)

sr = 20000;
dsr = 400;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
for o = 1:sf
    fname = files_short{o}     % filename
    ff = [inpdir fname];       % load
    load(ff)
    eeg = EEG.values';
    unit = unit.values;
    
    nqf = dsr / 2;      % filtering and downsampling unit
    flt = fir1(1024,40/nqf,'low');      % lowpass filtering at 40 Hz
    unit = filtfilt(flt,1,unit(1:const:end));
    
    eeg = eeg(1:const:end);     % downsampling EEG
    len = length(eeg);
    
    [pow_eeg,phase_eeg,f] = eegwavelet(unit,dsr);        % MPO wavelet
    
    H = zeros(size(pow_eeg,1),size(pow_eeg,1));         % cross-modulation
    MI = zeros(size(pow_eeg,1),size(pow_eeg,1));
    pbin = (phase_eeg * 180 / pi + 180) / 20;
    bin = fix(pbin) + 1;
    N = length(cnts);
    Hmax = log(N);
    for k1 = 1:size(pow_eeg,1)
        disp(k1)
        for k2 = 1:size(pow_eeg,1)
            A = zeros(1,N);
            for t = 1:N
                pw = pow_eeg(k2,:);
                A(t) = mean(pw(bin(k1,:)==t));
            end
            pb = A / sum(A);
            H(k1,k2) = -sum(pb.*log(pb));
            MI(k1,k2) = (Hmax - H(k1,k2)) / Hmax;
        end
    end
    H = figure;     % plot cross-modulation
    imagesc(fliplr(MI'))
    b_rescaleaxis('X',f(end:-1:1))
    b_rescaleaxis('Y',f)
    xlabel('Phase frequency [Hz]')
    ylabel('Amplitude frequency [Hz]')
    setappdata(gca,'scalex',f(end:-1:1))
    setappdata(gca,'scaley',f)
    b_zoomset_for_wavelet
    tt = fname(1:end-4);
    tt(tt=='_') = ' ';
    tt = [tt ' Unit crossmod'];
    title(tt)
    fn = [resdir tt '_UXMOD.fig'];       % save unit cross-modulation
    saveas(H,fn)
    fn = [resdir tt '_UXMOD.tiff'];
    saveas(H,fn)
        
    fn = [resdir fname(1:end-4) '_UXMOD.mat'];
    save(fn,'MI')

    waitbar(o/sf)
    close all
end



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
vrs = version;
if isequal(vrs(1:5),'7.4.0')
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
function [pow,phase,f,spf] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
dt = 1 / sr;
pad = 1;
spf = 5;
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
s = 1 ./ (50:-(1/spf):0.2) ./ c;
fperiod = c .* s;
f = 1 ./ fperiod;
param = -1;
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = wavelet(dat,dt,pad,s,mother,param);
pow = abs(wave).^2;
phase = angle(wave);

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
function b_zoomset_for_wavelet
%ZOOMSET_FOR_WAVELET    Assigns ZOOM_FOR_WAVELET as figure ButtonDownFcn.
%
%   See also ZOOM2 and ZOOMSET.

% Get handles
ax = gca;
plt = findobj(ax,'Type','line');
ptch = findobj(ax,'Type','patch');
im = findobj(ax,'Type','image');

% Set application data
if ~isempty(plt)
    x_data = get(plt,'XData');
    y_data = get(plt,'YData');
    setappdata(ax,'x_data',x_data)
    setappdata(ax,'y_data',y_data)
end
x_lim = get(ax,'XLim');
y_lim = get(ax,'YLim');
setappdata(ax,'x_lim',x_lim)
setappdata(ax,'y_lim',y_lim)

% Set ButtonDownFcn
set(ax,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
if ~isempty(plt)
    set(plt,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
elseif ~isempty(ptch)
    set(ptch,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
elseif ~isempty(im)
    set(im,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
end