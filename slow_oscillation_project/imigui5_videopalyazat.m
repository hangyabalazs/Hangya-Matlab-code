function varargout = imigui5_videopalyazat(varargin)
% IMIGUI5_VIDEOPALYAZAT M-file for imigui5_videopalyazat.fig
%      IMIGUI5_VIDEOPALYAZAT, by itself, creates a new IMIGUI5_VIDEOPALYAZAT or raises the existing
%      singleton*.
%
%      H = IMIGUI5_VIDEOPALYAZAT returns the handle to a new IMIGUI5_VIDEOPALYAZAT or the handle to
%      the existing singleton*.
%
%      IMIGUI5_VIDEOPALYAZAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMIGUI5_VIDEOPALYAZAT.M with the given input arguments.
%
%      IMIGUI5_VIDEOPALYAZAT('Property','Value',...) creates a new IMIGUI5_VIDEOPALYAZAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imigui5_videopalyazat_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imigui5_videopalyazat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imigui5_videopalyazat

% Last Modified by GUIDE v2.5 22-May-2008 16:09:33

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
pt = '40';
eg = '125';
ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\MIrawmaps_EEG' eg '_filt01_40_overlap10.mat'];
load(ff);nW=rI;
handles.nW = nW;
handles.Wmax = max(nW(:));
handles.Wmin = min(nW(:));
handles.WCOHrate = 1;     % 1 WCoh value / sec.

% Significance levels
ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\siglev_EEG' eg '.mat'];
% ff = [DATAPATH 'Ulbert\MImap\MIshiftcontrolfine.mat'];
load(ff)
handles.siglev = sl;

% Load MI shift data
ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\MIshiftfine.mat'];
% ff = [DATAPATH 'Ulbert\MImap\MIshiftcontrolfine.mat'];
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

% Original data
global DATADIR
fn = [DATADIR 'human_SO\oiti' pt '_mojzsis\grid\mat_EEG_' eg '\EEG_' eg '_0_65_rs.mat'];
% fn = 'F:\raw_data\human_SO\oiti37_lukacs\grid\mat\arteficial_control.mat';
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
set(handles.patch_handle,'CDataMapping','direct','EdgeColor','none')
axes(handles.axes3)
xlim([1 grid_horz])
ylim([1 grid_vert])
for k = 1:handles.chno
    [inx1 inx2] = gridind2sub(k,grid_horz);
    text(inx2,5-inx1,['{\it' num2str(k) '}'])
end
axis off
% axis([1 grid_horz 1 grid_vert])
% set(gca,'XLimMode','manual')
% set(gca,'YLimMode','manual')
% handles.patch2_handle = patch(X,Y,C);
% set(handles.patch2_handle,'CDataMapping','direct','EdgeColor','none')
handles.grid_horz = grid_horz;
handles.grid_vert = grid_vert;
handles.linecolormap = colormap('lines');
handles.colormap = colormap('jet');
handles.sizeC = size(C,2);

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
set(handles.figure1,'Units',fig_oldunits)
set(handles.axes3,'Units',ax_oldunits)

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
mx = size(handles.rI,3);
slst = 1 / (mx - 1);
set(handles.slider1,'SliderStep',[slst/3 slst])
% Get slider value
v = get(handles.slider1,'Value');

% Get time point
trI = handles.rI(:,:,floor(v)) + ...
    (handles.rI(:,:,ceil(v))-handles.rI(:,:,floor(v)))*(v-floor(v));

% profile off
% profile on -detail builtin -timer real
% tic

% Draw
grid_horz = handles.grid_horz;
grid_vert = handles.grid_vert;
numel_grid = grid_horz * grid_vert;
cm = handles.colormap;
sC = size(cm,1);
C = zeros(1,handles.sizeC);
next = 1;
for k = 1:numel_grid
    if k < grid_horz        % horizontal
        ci = rIclr(handles,trI(k,k+1),sC);
        C(1,next) = ci;
        next = next + 1;
    elseif k > grid_horz && rem(k,grid_horz) ~= 0 && k < numel_grid - grid_horz
        ci = rIclr(handles,trI(k,k+1),sC);
        C(1,next) = ci;
        next = next + 1;
    elseif k > numel_grid - grid_horz && k < numel_grid
        ci = rIclr(handles,trI(k,k+1),sC);
        C(1,next) = ci;
        next = next + 1;
    end
    if rem(k,grid_horz) == 1 && k < numel_grid - grid_horz          % vertical
        ci = rIclr(handles,trI(k,k+5),sC);
        C(1,next) = ci;
        next = next + 1;
    elseif rem(k,grid_horz) > 1 && rem(k,grid_horz) < grid_horz && k < numel_grid - grid_horz
        ci = rIclr(handles,trI(k,k+5),sC);
        C(1,next) = ci;
        next = next + 1;
    elseif rem(k,grid_horz) == 0 && k <= numel_grid - grid_horz
        ci = rIclr(handles,trI(k,k+5),sC);
        C(1,next) = ci;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 0 && k <= numel_grid - grid_horz     % diagonal
        ci = rIclr(handles,trI(k,k+6),sC);
        C(1,next) = ci;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 1 && k > grid_horz
        ci = rIclr(handles,trI(k,k-6),sC);
        C(1,next) = ci;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 1 && k <= numel_grid - grid_horz     % antidiagonal
        ci = rIclr(handles,trI(k,k+4),sC);
        C(1,next) = ci;
        next = next + 1;
    end
    if rem(k,grid_horz) ~= 0 && k > grid_horz
        ci = rIclr(handles,trI(k,k-4),sC);
        C(1,next) = ci;
        next = next + 1;
    end
end
set(handles.patch_handle,'CData',C)
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

% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = 5;
end