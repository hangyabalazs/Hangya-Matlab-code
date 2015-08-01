function hixmodrun
%HIXMODRUN   Runs cross-modulation on a sequence of files.
%   HIXMODRUN calculates and saves crossmodulation plots (see Tort et al.,
%   2008, PNAS).
%
%   See also WAVELET_NEW3.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = [DATAPATH 'Hajni\EEGMPO\mat\'];
resdir = [DATAPATH 'Hajni\EEGMPO_new\crossmod\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running HIXMODRUN...','Position',[360 250 275 50]);
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
    
    [pow_eeg,phase_eeg,f] = eegwavelet(eeg,dsr);        % EEG wavelet
%     [pow_unit,phase_unit,f] = eegwavelet(unit,dsr);        % unit wavelet
    
%     H = zeros(size(pow_eeg,1),size(pow_eeg,1));
%     MI = zeros(size(pow_eeg,1),size(pow_eeg,1));
%     for k1 = 1:size(pow_eeg,1)
%         for k2 = 1:size(pow_eeg,1)
%             [nm,bin] = histc(phase_eeg(k1,:)*180/pi,edges);   % phase histogram
%             nm = nm(1:end-1);
%             N = length(nm);
%             A = zeros(1,N);
%             for t = 1:N
%                 pw = pow_eeg(k2,:);
%                 A(t) = mean(pw(bin==t));
%             end
%             pb = A / sum(A);
%             H(k1,k2) = -sum(pb.*log(pb));
%             Hmax = log(N);
%             MI(k1,k2) = (Hmax - H(k1,k2)) / Hmax;
%         end
%     end
%     figure
%     imagesc(MI)
    
    H = zeros(size(pow_eeg,1),size(pow_eeg,1));
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
    figure
    imagesc(fliplr(MI'))
    b_rescaleaxis('X',f(end:-1:1))
    b_rescaleaxis('Y',f)
    setappdata(gca,'scalex',f(end:-1:1))
    setappdata(gca,'scaley',f)
    b_zoomset_for_wavelet
    fn = [resdir tt '_crossmod.fig'];       % save
    saveas(H,fn)
    
    figure
    imagesc(50:-0.2:0.2,0.2:0.2:50,MI')
    
%     H = zeros(size(pow_eeg,1),size(pow_eeg,1));
%     MI = zeros(size(pow_eeg,1),size(pow_eeg,1));
%     pbin = (phase_eeg * 180 / pi + 180) / 20;
%     bin = fix(pbin) + 1;
%     bw = (max(pow_eeg(:)) - min(pow_eeg(:))) / 100;
%     ppwbin = (pow_eeg - min(pow_eeg(:))) / bw;
%     pwbin = fix(ppwbin) + 1;
%     bin(bin>18) = 18;
%     pwbin(pwbin>100) = 100;
%     hpp = accumarray([bin(:) pwbin(:)],1,[18 100]);
%     hpp = hpp / sum(hpp(:));
%     N = length(cnts);
%     Hmax = log(N);
%     for k1 = 1:size(pow_eeg,1)
%         disp(k1)
%         for k2 = 1:size(pow_eeg,1)
%             A = zeros(1,N);
%             for t = 1:N
%                 pw = pow_eeg(k2,:);
%                 A(t) = mean(pw(bin(k1,:)==t));
%             end
%             pb = A / sum(A);
%             H(k1,k2) = -sum(pb.*log(pb));
%             MI(k1,k2) = (Hmax - H(k1,k2)) / Hmax;
%         end
%     end
%     figure
%     imagesc(MI)




    wave_cross = wave_eeg .* conj(wave_unit);       % crosswavelet
    im = abs(wave_cross) .^ 2;
    dt = 1 / dsr;
    wavetime = [0:len-1] * dt;
    
    pwind1_delta = find(f>4,1,'last');    % delta frequency band bounderies
    pwind2_delta = find(f<0.5,1,'first') - 1;
    pwind1_spindle = find(f>20,1,'last');    % spindle frequency band bounderies
    pwind2_spindle = find(f<7,1,'first');
    pwind1_gamma = find(f>40,1,'last');    % theta frequency band bounderies
    pwind2_gamma = find(f<30,1,'first');
    MeanDeltaPower = mean(im(pwind1_delta:pwind2_delta,:));
    MeanSpindlePower = mean(im(pwind1_spindle:pwind2_spindle,:));
    MeanGammaPower = mean(im(pwind1_gamma:pwind2_gamma,:));
    MaxDeltaPower = max(im(pwind1_delta:pwind2_delta,:));
    MaxSpindlePower = max(im(pwind1_spindle:pwind2_spindle,:));
    MaxGammaPower = max(im(pwind1_gamma:pwind2_gamma,:));
    CrossWaveMax = max(im);
    MeanDeltaPowerEEG = mean(im_eeg(pwind1_delta:pwind2_delta,:));
    MeanSpindlePowerEEG = mean(im_eeg(pwind1_spindle:pwind2_spindle,:));
    MeanGammaPowerEEG = mean(im_eeg(pwind1_gamma:pwind2_gamma,:));
    MaxDeltaPowerEEG = max(im_eeg(pwind1_delta:pwind2_delta,:));
    MaxSpindlePowerEEG = max(im_eeg(pwind1_spindle:pwind2_spindle,:));
    MaxGammaPowerEEG = max(im_eeg(pwind1_gamma:pwind2_gamma,:));
    MeanDeltaPowerUnit = mean(im_unit(pwind1_delta:pwind2_delta,:));
    MeanSpindlePowerUnit = mean(im_unit(pwind1_spindle:pwind2_spindle,:));
    MeanGammaPowerUnit = mean(im_unit(pwind1_gamma:pwind2_gamma,:));
    MaxDeltaPowerUnit = max(im_unit(pwind1_delta:pwind2_delta,:));
    MaxSpindlePowerUnit = max(im_unit(pwind1_spindle:pwind2_spindle,:));
    MaxGammaPowerUnit = max(im_unit(pwind1_gamma:pwind2_gamma,:));
        
    H = figure;         % plot crosswavelet
%     S1 = subplot(3,1,1);
    imagesc(10*log10(im))
    set(gca,'CLim',[5 40])
    mx = max(log(im(:)));
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
    ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
    text(xcoord,ycoord,num2str(mx),'Color',[1 1 1]);
    ff = round(f*100) / 100;
    time = round(wavetime*100) / 100;
    b_rescaleaxis('Y',ff)
    b_rescaleaxis('X',time)
    tt = fname(1:end-4);
    tt(tt=='_') = ' ';
    tt = [tt ' Crosswavelet'];
    title(tt)
    setappdata(gca,'scalex',wavetime)
    setappdata(gca,'scaley',f(f>0.5))
    b_zoomset_for_wavelet
%     S2 = subplot(3,1,2);
%     dneeg = (eeg - mean(eeg)) / std(eeg);
%     etime = 1:size(im,2);
%     plot(etime,dneeg)
%     xlim([etime(1) etime(end)])
%     hold on
%     flt = fir1(1024,[7 20]/nqf,'band');      % filtering in spindle range
%     fdneeg = filtfilt(flt,1,dneeg);
%     plot(etime,fdneeg,'r')
%     S3 = subplot(3,1,3);
%     dnunit = (unit - mean(unit)) / std(unit);
%     plot(etime,dnunit)
%     xlim([etime(1) etime(end)])
%     hold on
%     flt = fir1(1024,[7 20]/nqf,'band');      % filtering in spindle range
%     fdnunit = filtfilt(flt,1,dnunit);
%     plot(etime,fdnunit,'r')
%     setappdata(S1,'subplothandle',[S2 S3])
    fn = [resdir tt '_spindle.tiff'];       % save
    saveas(H,fn)
    fn = [resdir tt '_spindle.fig'];
    saveas(H,fn)
    H = figure;
    loglog(f(f>0.5),sum(im'))
    fn = [resdir tt '_xspect.fig'];       % save
    saveas(H,fn)
    
    H = figure;         % plot EEG wavelet
    imagesc(10*log10(im_eeg))
    set(gca,'CLim',[5 40])
    mx = max(log(im_eeg(:)));
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
    ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
    text(xcoord,ycoord,num2str(mx),'Color',[1 1 1]);
    ff = round(f*100) / 100;
    time = round(wavetime*100) / 100;
    b_rescaleaxis('Y',ff)
    b_rescaleaxis('X',time)
    tt = fname(1:end-4);
    tt(tt=='_') = ' ';
    tt = [tt ' EEG wavelet'];
    title(tt)
    setappdata(gca,'scalex',wavetime)
    setappdata(gca,'scaley',f)
    b_zoomset_for_wavelet
    fn = [resdir tt '_spindle.tiff'];       % save
    saveas(H,fn)
    H = figure;
    loglog(f(f>0.5),sum(im_eeg'))
    fn = [resdir tt '_eegspect.fig'];       % save
    saveas(H,fn)
    
    H = figure;         % plot unit wavelet
    imagesc(10*log10(im_unit))
    set(gca,'CLim',[5 40])
    mx = max(log(im_unit(:)));
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
    ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
    text(xcoord,ycoord,num2str(mx),'Color',[1 1 1]);
    ff = round(f*100) / 100;
    time = round(wavetime*100) / 100;
    b_rescaleaxis('Y',ff)
    b_rescaleaxis('X',time)
    tt = fname(1:end-4);
    tt(tt=='_') = ' ';
    tt = [tt ' Unit wavelet'];
    title(tt)
    setappdata(gca,'scalex',wavetime)
    setappdata(gca,'scaley',f)
    b_zoomset_for_wavelet
    fn = [resdir tt '_spindle.tiff'];       % save
    saveas(H,fn)
    H = figure;
    loglog(f(f>0.5),sum(im_unit'))
    fn = [resdir tt '_unitspect.fig'];       % save
    saveas(H,fn)
    
    fn = [resdir fname(1:end-4) '_CWVARS.mat'];
    save(fn,'MeanDeltaPower','MeanSpindlePower','MeanGammaPower',...
        'MaxDeltaPower','MaxSpindlePower','MaxGammaPower','CrossWaveMax',...
        'MeanDeltaPowerEEG','MeanSpindlePowerEEG','MeanGammaPowerEEG',...
        'MaxDeltaPowerEEG','MaxSpindlePowerEEG','MaxGammaPowerEEG',...
        'MeanDeltaPowerUnit','MeanSpindlePowerUnit','MeanGammaPowerUnit',...
        'MaxDeltaPowerUnit','MaxSpindlePowerUnit','MaxGammaPowerUnit')

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