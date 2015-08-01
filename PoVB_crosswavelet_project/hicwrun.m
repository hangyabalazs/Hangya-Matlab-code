function hicwrun
%HICWRUN   Runs crosswavelet on a sequence of files.
%   HICWRUN calculates and saves wavelet and crosswavelet images.
%
%   See also WAVELET_NEW3.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = [DATAPATH 'Hajni\EEGMPO\mat\'];
resdir = [DATAPATH 'Hajni\EEGMPO_new\crosswavelet\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running HICWRUN...','Position',[360 250 275 50]);
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
for o = 1:sf
    fname = files_short{o}     % filename
    ff = [inpdir fname];       % load
    load(ff)
    eeg = EEG.values';
    unit = unit.values;
    
    nqf = dsr / 2;      % filtering and downsampling unit
    flt = fir1(1024,40/nqf,'low');      % lowpass filtering on 100 Hz
    unit = filtfilt(flt,1,unit(1:const:end));
    
    eeg = eeg(1:const:end);     % downsampling EEG
    len = length(eeg);
    
    [wave_eeg,f] = eegwavelet(eeg,dsr);        % EEG wavelet
    im_eeg = abs(wave_eeg) .^ 2;
    [wave_unit,f] = eegwavelet(unit,dsr);        % unit wavelet
    im_unit = abs(wave_unit) .^ 2;
    
    wave_cross = wave_eeg .* conj(wave_unit);       % crosswavelet
    im = abs(wave_cross) .^ 2;
    dt = 1 / dsr;
    wavetime = [0:len-1] * dt;
    CrossWaveMax = max(im);
    
    H = figure;         % plot
    imagesc(im)
    mx = max(im(:));
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
    setappdata(gca,'scaley',f)
    b_zoomset_for_wavelet
    fn = [resdir tt '.tiff'];       % save
%     saveas(H,fn)
    
    H = figure;         % plot
    imagesc(im_eeg)
    mx = max(im_eeg(:));
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
    fn = [resdir tt '.tiff'];       % save
%     saveas(H,fn)
    
    H = figure;         % plot
    imagesc(im_unit)
    mx = max(im_unit(:));
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
    fn = [resdir tt '.tiff'];       % save
%     saveas(H,fn)
    
    fn = [resdir fname(1:end-4) '_CXMAX.mat'];
%     save(fn,'CrossWaveMax')

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
function [wave,f] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 1;
dj = 0.02;    
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