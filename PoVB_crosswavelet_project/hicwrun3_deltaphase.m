function hicwrun3_deltaphase
%HICWRUN3_DELTAPHASE   Runs crosswavelet on a sequence of files.
%   HICWRUN3_DELTAPHASE calculates delta-phase (0.5-4 Hz) dependence of
%   crosswavelet power in spindle range (7-20 Hz).
%
%   See also HICWRUN3.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = [DATAPATH 'Hajni\EEGMPO\mat2\'];
resdir = [DATAPATH 'Hajni\EEGMPO_nonorm\phasedep\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running HICWRUN3 DELTAPHASE...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
main(inpdir,resdir,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir,resdir,files_short,sf,wb)

cd(resdir)
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
    
    edges = -pi:2*pi/18:pi;     % phase histogram bin limits
    [mn_phase mn_cw mn mvl inx R p H Hc] = ...
        phasedep2(eeg',edges,MaxSpindlePower,dsr,0.5,4,flo,2,1);    % phase dependence of spindle power
    tt = fname(1:end-4);    % title & axis labels
    tt2 = tt;
    tt(tt=='_') = ' ';
    title(tt)
    xlabel('delta phase')
    ylabel('spindle power')
    fn = [resdir tt2 '_phasedep.fig'];       % save figures
    saveas(H,fn)
    fn = [resdir tt2 '_phasedep_ctrl.fig'];
    saveas(Hc,fn)
    
    xlsname = [resdir 'deltaphase_vs_MaxSpindlePower.xls'];   % write results to excel file
    if b_isfilename(xlsname)    % write Excel file
        [ntz mtz atz] = xlsread(xlsname,'sheet1');
        pref = size(atz,1) + 1;
    else
        pref = 1;
    end
    xlswrite(xlsname,{fname},'sheet1',['A' num2str(pref)])
    xlswrite(xlsname,mn*180/pi,'sheet1',['B' num2str(pref)])
    xlswrite(xlsname,mvl,'sheet1',['C' num2str(pref)])
    xlswrite(xlsname,R,'sheet1',['D' num2str(pref)])
    xlswrite(xlsname,p,'sheet1',['E' num2str(pref)])
    
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