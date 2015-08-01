function hicwrun3
%HICWRUN3   Runs crosswavelet on a sequence of files.
%   HICWRUN3 calculates and saves wavelet and crosswavelet images.
%
%   See also WAVELET_NEW3.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = [DATAPATH 'Hajni\EEGMPO\mat\'];
resdir = [DATAPATH 'Hajni\EEGMPO_nonorm\crosswavelet_spindle\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running HICWRUN3...','Position',[360 250 275 50]);
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
    try
        eeg = EEG.values';
    catch
        eeg = EEG1.values';
    end
    unit = unit.values;
    
    nqf = dsr / 2;      % filtering and downsampling unit
    if length(unit) / const > 3 * 1024
        flo = 1024;
    else
        flo = 512;
        disp('Filter order: 512')
    end
    flt = fir1(flo,40/nqf,'low');      % lowpass filtering at 40 Hz
    unit = filtfilt(flt,1,unit(1:const:end));
    unit = unit - mean(unit);
    eeg = eeg(1:const:end);     % downsampling EEG
    eeg = eeg - mean(eeg);
    len = length(eeg);
    sze = length(eeg);
    szu = length(unit);
    if ~isequal(sze,szu)
        szz = min(sze,szu);
        eeg = eeg(1:szz);
        unit = unit(1:szz);
    end
    
    [wave_eeg,f] = eegwavelet(eeg,dsr);        % EEG wavelet
    im_eeg = abs(wave_eeg) .^ 2;
    [wave_unit,f] = eegwavelet(unit,dsr);        % unit wavelet
    im_unit = abs(wave_unit) .^ 2;
    
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
    S1 = subplot(3,1,1);
    imagesc(im(1:pwind2_spindle,:))
    mx = max(max(im(1:pwind2_spindle,:)));
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
    setappdata(S1,'scalex',wavetime)
    setappdata(S1,'scaley',f(f>0.5))
    b_zoomset_for_wavelet
    S2 = subplot(3,1,2);
    dneeg = (eeg - mean(eeg)) / std(eeg);
    etime = 1:size(im,2);
    plot(etime,dneeg)
    xlim([etime(1) etime(end)])
    hold on
    flt = fir1(flo,[7 20]/nqf,'band');      % filtering in spindle range
    fdneeg = filtfilt(flt,1,dneeg);
    plot(etime,fdneeg,'r')
    S3 = subplot(3,1,3);
    dnunit = (unit - mean(unit)) / std(unit);
    plot(etime,dnunit)
    xlim([etime(1) etime(end)])
    hold on
    flt = fir1(flo,[7 20]/nqf,'band');      % filtering in spindle range
    fdnunit = filtfilt(flt,1,dnunit);
    plot(etime,fdnunit,'r')
    setappdata(S1,'subplothandle',[S2 S3])
    fn = [resdir tt '_spindle.tiff'];       % save
    saveas(H,fn)
    fn = [resdir tt '_spindle.fig'];
    saveas(H,fn)
    H = figure;
    loglog(f(f>0.5),sum(im'))
    fn = [resdir tt '_xspect.fig'];       % save
    saveas(H,fn)
    
    H = figure;         % plot EEG wavelet
    imagesc(im_eeg(1:pwind2_spindle,:))
    mx = max(max(im_eeg(1:pwind2_spindle,:)));
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
    imagesc(im_unit(1:pwind2_spindle,:))
    mx = max(max(im_unit(1:pwind2_spindle,:)));
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