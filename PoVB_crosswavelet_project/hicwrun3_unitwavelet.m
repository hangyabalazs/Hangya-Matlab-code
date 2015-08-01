function hicwrun3_unitwavelet
%HICWRUN3_UNITWAVELET   Runs wavelet on a sequence of files.
%   HICWRUN3_UNITWAVELET calculates and saves unit wavelet.
%
%   See also WAVELET_NEW3.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = [DATAPATH 'Hajni\EEGMPO\mat\'];
resdir = [DATAPATH 'Hajni\EEGMPO_new\unitwavelet\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running HICWRUN3 UNITWAVELET...','Position',[360 250 275 50]);
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
nqf = dsr / 2;
for o = 1:sf
    fname = files_short{o}     % filename
    ff = [inpdir fname];       % load
    load(ff)
    unit = unit.values;
    unit = unit(1:const:end);
    len = length(unit);
    
    [wave_unit,f] = eegwavelet(unit,dsr);        % unit wavelet
    im_unit = abs(wave_unit) .^ 2;
    dt = 1 / dsr;
    wavetime = (0:len-1) * dt;
    
    pwind1_highgamma = find(f>100,1,'last');    % high gamma frequency band bounderies
    pwind2_highgamma = find(f<50,1,'first');
    MeanHighGammaPowerUnit = mean(im_unit(pwind1_highgamma:pwind2_highgamma,:)); %#ok<NASGU>
    MaxHighGammaPowerUnit = max(im_unit(pwind1_highgamma:pwind2_highgamma,:)); %#ok<NASGU>
        
    H = figure;         % plot crosswavelet
    S1 = subplot(2,1,1);
    imagesc(im_unit(1:pwind2_highgamma,:))
    mx = max(im_unit(1:pwind2_highgamma,:));
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
    S2 = subplot(2,1,2);
    etime = 1:size(im_unit,2);
    dnunit = (unit - mean(unit)) / std(unit);
    plot(etime,dnunit)
    xlim([etime(1) etime(end)])
    hold on
    flt = fir1(1024,[50 100]/nqf,'band');      % filtering in high gamma range
    fdnunit = filtfilt(flt,1,dnunit);
    plot(etime,fdnunit,'r')
    setappdata(S1,'subplothandle',S2)
    fn = [resdir tt '_highgamma.tiff'];       % save
    saveas(H,fn)
    fn = [resdir tt '_highgamma.fig'];
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
    fn = [resdir tt '_log_highgamma.tiff'];       % save
    saveas(H,fn)
    H = figure;
    loglog(f(f>0.5),sum(im_unit'))
    fn = [resdir tt '_unitspect.fig'];       % save
    saveas(H,fn)
    
    fn = [resdir fname(1:end-4) '_CWVARS.mat'];
    save(fn,'MeanHighGammaPowerUnit','MaxHighGammaPowerUnit')

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