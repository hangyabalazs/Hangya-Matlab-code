function varargout = imigui2(varargin)
% IMIGUI2 M-file for imigui2.fig
%      IMIGUI2, by itself, creates a new IMIGUI2 or raises the existing
%      singleton*.
%
%      H = IMIGUI2 returns the handle to a new IMIGUI2 or the handle to
%      the existing singleton*.
%
%      IMIGUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMIGUI2.M with the given input arguments.
%
%      IMIGUI2('Property','Value',...) creates a new IMIGUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imigui2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imigui2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imigui2

% Last Modified by GUIDE v2.5 22-May-2008 16:09:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imigui2_OpeningFcn, ...
                   'gui_OutputFcn',  @imigui2_OutputFcn, ...
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
function imigui2_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for imigui2
handles.output = hObject;

% Load mutual information map
global DATAPATH
ff = [DATAPATH 'Ulbert\MImap\MIrawmaps_filtered40.mat'];
load(ff)
handles.rI = rI;
handles.Imax = max(rI(:));
handles.Imin = min(rI(:));
handles.MIrate = 1;     % 1 MI value / sec.

% Load wavelet coherence map
global DATAPATH
ff = [DATAPATH 'Ulbert\Wcoh\WCOHmaps.mat'];
load(ff);
handles.nW = nW;
handles.Wmax = max(nW(:));
handles.Wmin = min(nW(:));
handles.WCOHrate = 1;     % 1 WCoh value / sec.

% Original data
fn = 'F:\raw_data\human_SO\oiti37_lukacs\grid\mat\EEG12_4x5grid_18i_rs.mat';
load(fn)
handles.data = data;
handles.normdata = data - repmat(mean(data),size(data,1),1);
handles.normmax = max(handles.normdata(:));
handles.normmin = min(handles.normdata(:));
[pthname fname ext] = fileparts(fn);
handles.datname = [pthname '\' fname '.dat'];
handles.nfoname = [pthname '\' fname '.nfo'];
fid = fopen(handles.nfoname);
handles.chno = fread(fid,1,'double');
handles.datalen = fread(fid,1,'double');
handles.frate = 1000;       % data sampled on 1000 Hz

% Initialize slider
mx = size(handles.rI,3);
slst = 1 / (mx - 1);
set(handles.slider1,'Min',1,'Max',mx,'Value',1,'SliderStep',[slst/10 slst/5])

% Initialize grid
axes(handles.axes1)
axis off
grid_horz = 5;
grid_vert = 4;
axis([1 grid_horz 1 grid_vert])
set(gca,'XLimMode','manual')
set(gca,'YLimMode','manual')
numel_grid = grid_horz * grid_vert;
ep = zeros(grid_vert,grid_horz,2);      % electrode points
ep(:,:,1) = repmat((1:grid_horz),grid_vert,1);
ep(:,:,2) = repmat((grid_vert:-1:1)',1,grid_horz);
X = [];
Y = [];
C = [];
for k = 1:numel_grid
    [y x] = ind2sub([grid_horz,grid_vert],k);
    if k < grid_horz
        px = [ep(x,y,1) ep(x,y+1,1) ep(x,y,1)+1/2 ep(x,y,1) ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y+1,2) ep(x,y,2)-1/4 ep(x,y,2) ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0.25 0 0];
    elseif k > grid_horz && rem(k,grid_horz) ~= 0 && k < numel_grid - grid_horz
        px = [ep(x,y,1) ep(x,y,1)+1/2 ep(x,y+1,1) ep(x,y,1)+1/2 ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)+1/4 ep(x,y+1,2) ep(x,y,2)-1/4 ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0 0.5 0];
    elseif k > numel_grid - grid_horz && k < numel_grid
        px = [ep(x,y,1) ep(x,y+1,1) ep(x,y,1)+1/2 ep(x,y,1) ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y+1,2) ep(x,y,2)+1/4 ep(x,y,2) ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0 0 0.5];
    end
    if rem(k,grid_horz) == 1 && k < numel_grid - grid_horz
        px = [ep(x,y,1) ep(x+1,y,1) ep(x,y,1)+1/4 ep(x,y,1) ep(x,y,1)];
        py = [ep(x,y,2) ep(x+1,y,2) ep(x,y,2)-1/2 ep(x,y,2) ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0.5 0 0];
    elseif rem(k,grid_horz) > 1 && rem(k,grid_horz) < grid_horz && k < numel_grid - grid_horz
        px = [ep(x,y,1) ep(x,y,1)-1/4 ep(x+1,y,1) ep(x,y,1)+1/4 ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)-1/2 ep(x+1,y,2) ep(x,y,2)-1/2 ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0 0.25 0];
    elseif rem(k,grid_horz) == 0 && k <= numel_grid - grid_horz
        px = [ep(x,y,1) ep(x,y,1)-1/4 ep(x+1,y,1) ep(x,y,1) ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)-1/2 ep(x+1,y,2) ep(x,y,2) ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0 0 0.25];
    end
    if rem(k,grid_horz) ~= 0 && k <= numel_grid - grid_horz
        px = [ep(x,y,1) ep(x,y,1)+1/4 ep(x,y,1)+1/2 ep(x,y,1)+1/2 ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)-1/2 ep(x,y,2)-1/2  ep(x,y,2)-1/4 ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0.75 0 0];
    end
    if rem(k,grid_horz) ~= 1 && k > grid_horz
        px = [ep(x,y,1) ep(x,y,1)-1/4 ep(x,y,1)-1/2 ep(x,y,1)-1/2 ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)+1/2 ep(x,y,2)+1/2 ep(x,y,2)+1/4 ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0.75 0.75 0];
    end
    if rem(k,grid_horz) ~= 1 && k <= numel_grid - grid_horz
        px = [ep(x,y,1) ep(x,y,1)-1/4 ep(x,y,1)-1/2 ep(x,y,1)-1/2 ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)-1/2 ep(x,y,2)-1/2 ep(x,y,2)-1/4 ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0 0.75 0];
    end
    if rem(k,grid_horz) ~= 0 && k > grid_horz
        px = [ep(x,y,1) ep(x,y,1)+1/4 ep(x,y,1)+1/2 ep(x,y,1)+1/2 ep(x,y,1)];
        py = [ep(x,y,2) ep(x,y,2)+1/2 ep(x,y,2)+1/2 ep(x,y,2)+1/4 ep(x,y,2)];
        X = [X px'];
        Y = [Y py'];
        C(1,end+1,1:3) = [0 0 0.75];
    end
end
handles.patch_handle = patch(X,Y,C);
set(handles.patch_handle,'CDataMapping','direct')
axes(handles.axes3)
axis off
axis([1 grid_horz 1 grid_vert])
set(gca,'XLimMode','manual')
set(gca,'YLimMode','manual')
handles.patch2_handle = patch(X,Y,C);
set(handles.patch2_handle,'CDataMapping','direct')
handles.grid_horz = grid_horz;
handles.grid_vert = grid_vert;
handles.linecolormap = colormap('lines');
handles.colormap = colormap('jet');
handles.sizeC = size(C,2);

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
I = imagesc(rand(4,5)*1000,[handles.normmin handles.normmax]);
axis off
handles.amaphandle = I;

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------
function varargout = imigui2_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------
function slider1_Callback(hObject, eventdata, handles)

% Get slider value
v = get(handles.slider1,'Value');

% Get time point
trI = handles.rI(:,:,floor(v)) + ...
    (handles.rI(:,:,ceil(v))-handles.rI(:,:,floor(v)))*(v-floor(v));
tnW = handles.nW(:,:,floor(v)) + ...
    (handles.nW(:,:,ceil(v))-handles.nW(:,:,floor(v)))*(v-floor(v));

% Draw
% profile off
% profile on -detail builtin -timer real
% tic
grid_horz = handles.grid_horz;
grid_vert = handles.grid_vert;
numel_grid = grid_horz * grid_vert;
cm = handles.colormap;
sC = size(cm,1);
C = zeros(1,handles.sizeC);
Cw = zeros(1,handles.sizeC);
next = 1;
for k = 1:numel_grid
    if k < grid_horz        % horizontal
        ci = rIclr(handles,trI(k,k+1),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+1),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    elseif k > grid_horz && rem(k,grid_horz) ~= 0 && k < numel_grid - grid_horz
        ci = rIclr(handles,trI(k,k+1),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+1),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    elseif k > numel_grid - grid_horz && k < numel_grid
        ci = rIclr(handles,trI(k,k+1),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+1),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    end
    if rem(k,grid_horz) == 1 && k < numel_grid - grid_horz          % vertical
        ci = rIclr(handles,trI(k,k+5),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+5),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    elseif rem(k,grid_horz) > 1 && rem(k,grid_horz) < grid_horz && k < numel_grid - grid_horz
        ci = rIclr(handles,trI(k,k+5),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+5),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    elseif rem(k,grid_horz) == 0 && k <= numel_grid - grid_horz
        ci = rIclr(handles,trI(k,k+5),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+5),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 0 && k <= numel_grid - grid_horz     % diagonal
        ci = rIclr(handles,trI(k,k+6),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+6),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 1 && k > grid_horz
        ci = rIclr(handles,trI(k,k-6),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k-6),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 1 && k <= numel_grid - grid_horz     % antidiagonal
        ci = rIclr(handles,trI(k,k+4),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k+4),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 0 && k > grid_horz
        ci = rIclr(handles,trI(k,k-4),sC);
        C(1,next) = ci;
        ciw = nWclr(handles,tnW(k,k-4),sC);
        Cw(1,next) = ciw;
        next = next + 1;
    end
end
set(handles.patch_handle,'CData',C)
set(handles.patch2_handle,'CData',Cw)
% toc
% profile viewer

% Visualize original data
% profile off
% profile on -detail builtin -timer real
% tic
timep = v / handles.MIrate * handles.frate;
timep = round(timep);
time1 = timep - 6 * handles.frate;
time2 = timep + 5 * handles.frate;      % 11 sec. window
time3 = timep - 1 * handles.frate;      % 1 sec. window: time3-timep
time = linspace(time1,time2,11*handles.frate) / handles.frate;
axislim = [time1/handles.frate time2/handles.frate -2 22];
if time1 < 0
    pretag = zeros(1,-time1);
    posttag = [];
    time1 = 0;
elseif time2 > handles.datalen
    pretag = [];
    posttag = zeros(1,time2-handles.datalen);
    time2 = handles.datalen;
else
    pretag = [];
    posttag = [];
end
amp = zeros(1,handles.chno);
for k = 1:handles.chno
    chdata = handles.normdata(time1+1:time2,k);
    chdata = chdata / 700;     % adjust data
    chdata = [pretag chdata' posttag];
    set(handles.plot_handle(k),'XData',time,'YData',-chdata+k,'Color',...
        handles.linecolormap(k,:))
    chdata2 = handles.normdata(time3+1:timep,k);
    amp(k) = max(chdata2) - min(chdata2);
end
axis(handles.axes2,axislim)

% Amplitude map
set(handles.amaphandle,'CData',reshape(amp,handles.grid_horz,handles.grid_vert)')
% toc
% profile viewer

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