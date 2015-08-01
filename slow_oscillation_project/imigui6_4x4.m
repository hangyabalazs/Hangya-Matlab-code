function varargout = imigui6_4x4(varargin)
% IMIGUI6_4X4 M-file for imigui6_4x4.fig
%      IMIGUI6_4X4, by itself, creates a new IMIGUI6_4X4 or raises the existing
%      singleton*.
%
%      H = IMIGUI6_4X4 returns the handle to a new IMIGUI6_4X4 or the handle to
%      the existing singleton*.
%
%      IMIGUI6_4X4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMIGUI6_4X4.M with the given input arguments.
%
%      IMIGUI6_4X4('Property','Value',...) creates a new IMIGUI6_4X4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imigui5_4x4_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imigui6_4x4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imigui6_4x4

% Last Modified by GUIDE v2.5 28-Apr-2010 17:20:39

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
pt = 'N1';
eg = '157';
bi = [DATAPATH 'Ulbert\reko\reko3_OITI' pt '.jpg'];
ff = [DATAPATH 'Ulbert\OITI_' pt '_EEG_' eg '\MImap\MIrawmaps_EEG157_filt01_40_overlap10.mat'];
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
fn = [DATADIR 'human_SO\oiti' lower(pt) '_wittner\grid\mat_EEG_' eg '\EEG_' eg '_1900_1965_rs.mat'];
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
        arrh(x,y) = annotation('arrow',[annp(x,1) annp(y,1)],[annp(x,2) annp(y,2)],'Color',[0 0 0]);
        arrh(y,x) = annotation('arrow',[annp(y,1) annp(x,1)],[annp(y,2) annp(x,2)],'Color',[0 0 0]);
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
        if ~isnan(troIM(x,y)) && trIML(x,y) > 1 && troIM(x,y) > handles.siglev(4,2)   % sig. lev.: 0.0001
            set([handles.arrow_handle(x,y) handles.arrow_handle2(x,y)],'LineWidth',(1-trIML(x,y)/1001)*10+eps,...
                'Color',cm(round(trIM(x,y)*63)+1,:),'Visible','on');
        else
            set([handles.arrow_handle(x,y) handles.arrow_handle2(x,y)],'Visible','off')
        end
        if ~isnan(troIM(y,x)) && trIML(y,x) > 1 && troIM(y,x) > handles.siglev(4,2)
            set([handles.arrow_handle(y,x) handles.arrow_handle2(y,x)],'LineWidth',(1-trIML(y,x)/1001)*10+eps,...
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
    inx2 = nm_cols;
end

% -------------------------------------------------------------------------
function cG = gridloc_oiti31_virag(chnum)

cG = zeros(chnum,2);
cG(1,:) = [162 253];
cG(2,:) = [173 240];
cG(3,:) = [191 220];
cG(4,:) = [216 199];
cG(5,:) = [191 275];
cG(6,:) = [205 261];
cG(7,:) = [222 244];
cG(8,:) = [245 219];
cG(9,:) = [216 298];
cG(10,:) = [225 285];
cG(11,:) = [246 268];
cG(12,:) = [268 244];
cG(13,:) = [241 319];
cG(14,:) = [252 308];
cG(15,:) = [271 290];
cG(16,:) = [293 267];

% -------------------------------------------------------------------------
function cG = gridloc_oiti39_gaal(chnum)

cG = zeros(chnum,2);
cG(1,:) = [690 232];
cG(2,:) = [647 209];
cG(3,:) = [606 186];
cG(4,:) = [564 165];
cG(5,:) = [715 212];
cG(6,:) = [675 188];
cG(7,:) = [635 164];
cG(8,:) = [594 142];
cG(9,:) = [734 196];
cG(10,:) = [694 173];
cG(11,:) = [659 149];
cG(12,:) = [618 126];
cG(13,:) = [746 183];
cG(14,:) = [713 157];
cG(15,:) = [679 137];
cG(16,:) = [637 112];

% -------------------------------------------------------------------------
function cG = gridloc_oitin1_wittner(chnum)

cG = zeros(chnum,2);
cG(1,:) = [454 355];
cG(2,:) = [422 384];
cG(3,:) = [384 408];
cG(4,:) = [352 438];
cG(5,:) = [427 327];
cG(6,:) = [385 350];
cG(7,:) = [350 378];
cG(8,:) = [320 406];
cG(9,:) = [396 285];
cG(10,:) = [358 313];
cG(11,:) = [318 348];
cG(12,:) = [286 372];
cG(13,:) = [373 245];
cG(14,:) = [334 274];
cG(15,:) = [289 301];
cG(16,:) = [260 337];