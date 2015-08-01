function varargout = imigui(varargin)
% IMIGUI M-file for imigui.fig
%      IMIGUI, by itself, creates a new IMIGUI or raises the existing
%      singleton*.
%
%      H = IMIGUI returns the handle to a new IMIGUI or the handle to
%      the existing singleton*.
%
%      IMIGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMIGUI.M with the given input arguments.
%
%      IMIGUI('Property','Value',...) creates a new IMIGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imigui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imigui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imigui

% Last Modified by GUIDE v2.5 22-May-2008 16:09:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imigui_OpeningFcn, ...
                   'gui_OutputFcn',  @imigui_OutputFcn, ...
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
function imigui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for imigui
handles.output = hObject;

% Load mutual information map
global DATAPATH
ff = [DATAPATH 'Ulbert\MImap\MImaps.mat'];
load(ff)
handles.rI = rI;
handles.max = max(rI(:));
handles.min = min(rI(:));

% Initialize slider
mx = size(handles.rI,3);
slst = 1 / (mx - 1);
set(handles.slider1,'Min',1,'Max',mx,'Value',1,'SliderStep',[slst/10 slst/5])

% Axis off
axes(handles.axes1)
axis off

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------
function varargout = imigui_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------
function slider1_Callback(hObject, eventdata, handles)

% Get slider value
v = get(handles.slider1,'Value');

% Get time point
% mx = size(handles.rI,3);
% inx = min(max(round(v*mx),1),mx);
trI = handles.rI(:,:,floor(v)) + ...
    (handles.rI(:,:,ceil(v))-handles.rI(:,:,floor(v)))*(v-floor(v));

% Draw
profile off
profile on -detail builtin -timer real
tic
axes(handles.axes1)
% imagesc(imimap(trI),[1 2.1])
imagesc(imimap(trI),[handles.min handles.max])
axis off
toc
profile viewer

% -------------------------------------------------------------------------
function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


