function varargout = imigui6(varargin)
% IMIGUI6 M-file for imigui6.fig
%      IMIGUI6, by itself, creates a new IMIGUI6 or raises the existing
%      singleton*.
%
%      H = IMIGUI6 returns the handle to a new IMIGUI6 or the handle to
%      the existing singleton*.
%
%      IMIGUI6('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMIGUI6.M with the given input arguments.
%
%      IMIGUI6('Property','Value',...) creates a new IMIGUI6 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imigui5_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imigui6_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imigui6

% Last Modified by GUIDE v2.5 06-May-2010 11:25:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imigui3_OpeningFcn, ...
                   'gui_OutputFcn',  @imigui3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% -------------------------------------------------------------------------
function imigui3_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for imigui3
handles.output = hObject;

% Load wavelet coherence map
global DATAPATH
pt = '37';
eg = '12c';
bi = [DATAPATH 'Ulbert\reko\reko4_OITI' pt '.jpg'];
ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\MIrawmaps_EEG' eg '_filt01_40_overlap10.mat'];
load(ff);nW=rI;
handles.nW = nW;
handles.Wmax = max(nW(:));
handles.Wmin = min(nW(:));
handles.WCOHrate = 1;     % 1 WCoh value / sec.

% Significance levels
% ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\siglev_EEG' eg '.mat'];
ff = [DATAPATH 'Ulbert\EEG_12\MImap_SNRcontrol\siglev_EEG12.mat'];
load(ff)
handles.siglev = sl;

% Load MI shift data
% ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\MIshiftfine.mat'];
ff = [DATAPATH 'Ulbert\EEG_12\MImap_SNRcontrol\MIshiftfine.mat'];
load(ff)
handles.rIMaxLoc = rIMaxLoc;
handles.rIMax = rIMax;
handles.rIMaxmax = max(rIMax(:));
handles.rIMaxmin = min(rIMax(rIMax>handles.siglev(4,2)));
handles.normrIMax = (rIMax - handles.rIMaxmin) / (handles.rIMaxmax - handles.rIMaxmin);
handles.rIMaxrate = 10;     % 10 rIMax value / sec.

% Mutual information map (maximal MI!)
pr = permute(rIMax,[2 1 3]);
rpr(:,:,:,1) = pr;
rpr(:,:,:,2) = rIMax;
rmx = max(rpr,[],4);
handles.rI = rmx;
handles.Imax = max(rmx(:));
handles.Imin = min(rmx(:));
handles.MIrate = 10;     % 10 MI value / sec.

% Loading the MR image
grid_vert = 4;
if isequal(size(rIMax,1),20)
    grid_horz = 5;
else
    grid_horz = 4;
end
chnum = grid_horz * grid_vert;
cG = eval(['gridloc_oiti' pt '(chnum);']);    % electrode coordinates
I = imread(bi);         % read MR

% Original data
global DATADIR
% fn = [DATADIR 'human_SO\oiti' pt '_lukacs\grid\mat_EEG_' eg '\EEG_' eg '_4100_4160_rs.mat'];
fn = 'F:\raw_data\human_SO\oiti37_lukacs\grid\mat_EEG_12\arteficial_control_snr.mat';
load(fn)
data = data(1:end-7,:);
handles.data = data;
handles.normdata = data - repmat(mean(data),size(data,1),1);
handles.normmax = max(handles.normdata(:));
handles.normmin = min(handles.normdata(:));
[pthname fname ext] = fileparts(fn);
handles.datname = [pthname '\' fname '.dat'];
handles.nfoname = [pthname '\' fname '.nfo'];
handles.chno = size(handles.data,2);
handles.datalen = size(handles.data,1);
handles.frate = 1000;       % data sampled on 1000 Hz

% Initialize slider
mx = size(handles.rI,3);
slst = 1 / (mx - 1);
set(handles.slider1,'Min',1,'Max',mx,'Value',1,'SliderStep',[slst/10 slst])

% Initialize plot
axes(handles.axes1)
imagesc(I)
axis off

axes(handles.axes3)
xlim([1 grid_horz])
ylim([1 grid_vert])
for k = 1:handles.chno
    [inx1 inx2] = gridind2sub(k,grid_horz);
    text(inx2,5-inx1,['{\it' num2str(k) '}'])
end
axis off
handles.grid_horz = grid_horz;
handles.grid_vert = grid_vert;
handles.linecolormap = colormap('lines');
handles.colormap = colormap('jet');

% Annotation arrows
fig_oldunits = get(handles.figure1,'Units');      % rectangle annotation
set(handles.figure1,'Units','normalized')
ax_oldunits = get(handles.axes3,'Units');
set(handles.axes3,'Units','normalized')
ps = get(handles.axes3,'Position');
annp = zeros(handles.chno,2);      % annotaion arrow endpoints
for k = 1:handles.chno
    annp(k,1) = ps(1) + ps(3) - rem(handles.chno-k,grid_horz) * ps(3) / (grid_horz - 1);
    annp(k,2) = ps(2) + ps(4) - (ceil(k/grid_horz) - 1) * ps(3) / (grid_vert - 1);
end
arrh = zeros(handles.chno,handles.chno);      % arrow handles
for x = 1:handles.chno
    for y = x+1:handles.chno
        arrh(x,y) = annotation('arrow',[annp(x,1) annp(y,1)],[annp(x,2) annp(y,2)],'Color',[0 0 0]);
        arrh(y,x) = annotation('arrow',[annp(y,1) annp(x,1)],[annp(y,2) annp(x,2)],'Color',[0 0 0]);
    end
end
handles.arrow_handle = arrh;
set(handles.axes3,'Units',ax_oldunits)

ax_oldunits = get(handles.axes1,'Units');
set(handles.axes1,'Units','normalized')
ps = get(handles.axes1,'Position');
annp = zeros(handles.chno,2);      % annotaion arrow endpoints
for k = 1:handles.chno
    annp(k,1) = ps(1) + ps(3) / size(I,2) * cG(k,1);
    annp(k,2) = ps(2) + ps(4) - ps(4) / size(I,1) * cG(k,2);
end
arrh = zeros(handles.chno,handles.chno);      % arrow handles
for x = 1:handles.chno
    for y = x+1:handles.chno
        arrh(x,y) = annotation('arrow',[annp(x,1) annp(y,1)],[annp(x,2) annp(y,2)],'Color',[0 0 0],'HeadLength',5,'HeadWidth',5);
        arrh(y,x) = annotation('arrow',[annp(y,1) annp(x,1)],[annp(y,2) annp(x,2)],'Color',[0 0 0],'HeadLength',5,'HeadWidth',5);
    end
end
handles.arrow_handle2 = arrh;
set(handles.figure1,'Units',fig_oldunits)
set(handles.axes1,'Units',ax_oldunits)

% Raw data axes
axes(handles.axes2)
hold on
time1 = -5 * handles.frate;
time2 = 6 * handles.frate;      % 11 sec. window
time = linspace(time1,time2,11*handles.frate) / handles.frate;
axislim = [time1/handles.frate time2/handles.frate -2 22];
handles.plot_handle = zeros(1,handles.chno);
for k = 1:handles.chno
    handles.plot_handle(k) = plot(time,zeros(size(time))+k);
end
axis(axislim)

set(gca,'YDir','reverse')      % reverse y axis
fig_oldunits = get(handles.figure1,'Units');      % rectangle annotation
set(handles.figure1,'Units','normalized')
ax_oldunits = get(handles.axes2,'Units');
set(handles.axes2,'Units','normalized')
ps = get(handles.axes2,'Position');
lnx = ps(1) + 5 * ps(3) / 11;
wd = ps(3) / 11;
lny = ps(2) + 0.02;
hg = ps(4) - 0.04;
annotation('rectangle',[lnx lny wd hg]);
set(handles.figure1,'Units',fig_oldunits)
set(handles.axes2,'Units',ax_oldunits)

% Amplitudo map
axes(handles.axes4)
I = imagesc(rand(4,5)*1000,[0 handles.normmax-handles.normmin]);
axis off
handles.amaphandle = I;

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------
function varargout = imigui3_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------
function slider1_Callback(hObject, eventdata, handles)

% Get slider value
v = get(handles.slider1,'Value');

% Get time point
trIML = handles.rIMaxLoc(:,:,floor(v));
trIM = handles.normrIMax(:,:,floor(v));
troIM = handles.rIMax(:,:,floor(v));
cm = handles.colormap;

% Visualize original data
timep = 1000 + (v - 1) / handles.MIrate * handles.frate;
timep = round(timep);
time1 = timep - 6 * handles.frate;
time2 = timep + 5 * handles.frate;      % 11 sec. window
time3 = timep - 1 * handles.frate;      % 1 sec. window: time3-timep
time = linspace(time1,time2,11*handles.frate) / handles.frate;
axislim = [time1/handles.frate time2/handles.frate -2 22];
if time1 < 0
    pretag = zeros(1,-time1);
    time1 = 0;
else
    pretag = [];
end
if time2 > handles.datalen
    posttag = zeros(1,time2-handles.datalen);
    time2 = handles.datalen;
else
    posttag = [];
end
amp = zeros(1,handles.chno);
for k = 1:handles.chno
    chdata = handles.normdata(time1+1:time2,k);
    chdata = chdata / (handles.normmax - handles.normmin) * 2;     % adjust data
    chdata = [pretag chdata' posttag];
    set(handles.plot_handle(k),'XData',time,'YData',-chdata+k,'Color',...
        handles.linecolormap(k,:))
    chdata2 = handles.normdata(time3+1:timep,k);
    amp(k) = max(chdata2) - min(chdata2);
end
axis(handles.axes2,axislim)

% Amplitude map
set(handles.amaphandle,'CData',reshape(amp,handles.grid_horz,handles.grid_vert)')

% Plot MI-shift
for x = 1:handles.chno
    for y = x+1:handles.chno  % color: strength; width: speed
        if ~isnan(troIM(x,y)) && trIML(x,y) > 1 &&  troIM(x,y) > handles.siglev(4,2,x,y)   % sig. lev.: 0.0001
            set(handles.arrow_handle(x,y),'LineWidth',(1-trIML(x,y)/1001)*10+eps,...
                'Color',cm(round(trIM(x,y)*63)+1,:),'Visible','on');
            set(handles.arrow_handle2(x,y),'LineWidth',(1-trIML(x,y)/1001)*3+eps,...
                'Color',cm(round(trIM(x,y)*63)+1,:),'Visible','on');
        else
            set([handles.arrow_handle(x,y) handles.arrow_handle2(x,y)],'Visible','off')
        end
        if ~isnan(troIM(y,x)) && trIML(y,x) > 1 &&  troIM(y,x) > handles.siglev(4,2,y,x)
            set(handles.arrow_handle(y,x),'LineWidth',(1-trIML(y,x)/1001)*10+eps,...
                'Color',cm(round(trIM(y,x)*63)+1,:),'Visible','on');
            set(handles.arrow_handle2(y,x),'LineWidth',(1-trIML(y,x)/1001)*3+eps,...
                'Color',cm(round(trIM(y,x)*63)+1,:),'Visible','on');
        else
            set([handles.arrow_handle(y,x) handles.arrow_handle2(y,x)],'Visible','off')
        end
    end
end

% -------------------------------------------------------------------------
function c = rIclr(handles,t,sC)

pci = t;
pci2 = (pci - handles.Imin) / (handles.Imax - handles.Imin);   % 0 - 1
pci3 = pci2 * sC + 0.5;     % 0.5 - 64.5
pci4 = round(pci3);         % 1 - 64, 65 if pci3 == 64.5
c = min(pci4,sC);           % 1 - 64

% -------------------------------------------------------------------------
function c = nWclr(handles,t,sC)

pci = t;
pci2 = (pci - handles.Wmin) / (handles.Wmax - handles.Wmin);   % 0 - 1
pci3 = pci2 * sC + 0.5;     % 0.5 - 64.5
pci4 = round(pci3);         % 1 - 64, 65 if pci3 == 64.5
c = min(pci4,sC);           % 1 - 64

% -------------------------------------------------------------------------
function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% -------------------------------------------------------------------------
function y = var2(x,xbar)

n = length(x);
denom = n - 1;      % the unbiased estimator: divide by (n-1)
x0 = x - xbar;
y = sum(x0.^2) ./ denom;

% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = 5;
end

% -------------------------------------------------------------------------
function cG = gridloc_oiti40(chnum)

cG = zeros(chnum,2);
cG(1,:) = [361 194];
cG(2,:) = [384 236];
cG(3,:) = [403 274];
cG(4,:) = [425 314];
cG(5,:) = [444 347];
cG(6,:) = [313 218];
cG(7,:) = [333 249];
cG(8,:) = [355 285];
cG(9,:) = [375 328];
cG(10,:) = [397 366];
cG(11,:) = [273 237];
cG(12,:) = [294 269];
cG(13,:) = [318 304];
cG(14,:) = [337 344];
cG(15,:) = [361 379];
cG(16,:) = [227 245];
cG(17,:) = [250 284];
cG(18,:) = [273 321];
cG(19,:) = [296 361];
cG(20,:) = [319 401];

% -------------------------------------------------------------------------
function cG = gridloc_oiti37(chnum)

cG = zeros(chnum,2);
cG(1,:) = [784 260];
cG(2,:) = [772 280];
cG(3,:) = [754 313];
cG(4,:) = [732 349];
cG(5,:) = [706 379];
cG(6,:) = [758 227];
cG(7,:) = [744 251];
cG(8,:) = [724 281];
cG(9,:) = [705 316];
cG(10,:) = [673 344];
cG(11,:) = [733 199];
cG(12,:) = [716 223];
cG(13,:) = [689 248];
cG(14,:) = [669 287];
cG(15,:) = [638 314];
cG(16,:) = [705 173];
cG(17,:) = [688 195];
cG(18,:) = [665 223];
cG(19,:) = [636 256];
cG(20,:) = [606 283];


% -------------------------------------------------------------------------
function figure1_KeyPressFcn(hObject, eventdata, handles)

% Get application data
H = gcf;
inp = get(H,'CurrentCharacter');

% Play
handles = guidata(gcf);
switch inp
    case 'p'    % call slider
        for k = 1:700
            v = get(handles.slider1,'Value');
            set(handles.slider1,'Value',v+0.2);
            slider1_Callback(H,[],handles)
            pause(0.01)
        end
end
disp(num2str(v))
disp('That was it.')