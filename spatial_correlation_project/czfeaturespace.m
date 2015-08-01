function varargout = czfeaturespace(varargin)
% CZFEATURESPACE M-file for czfeaturespace.fig
%      CZFEATURESPACE, by itself, creates a new CZFEATURESPACE or raises the existing
%      singleton*.
%
%      H = CZFEATURESPACE returns the handle to a new CZFEATURESPACE or the handle to
%      the existing singleton*.
%
%      CZFEATURESPACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CZFEATURESPACE.M with the given input arguments.
%
%      CZFEATURESPACE('Property','Value',...) creates a new CZFEATURESPACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before czfeaturespace_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to czfeaturespace_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help czfeaturespace

% Last Modified by GUIDE v2.5 18-Feb-2010 10:22:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @czfeaturespace_OpeningFcn, ...
                   'gui_OutputFcn',  @czfeaturespace_OutputFcn, ...
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
function czfeaturespace_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for czfeaturespace
handles.output = hObject;

% get settings
global FEATURE_VECTOR
if ~isempty(FEATURE_VECTOR)
    if FEATURE_VECTOR.pc1
        set(handles.checkbox_pc1,'Value',1)
    else
        set(handles.checkbox_pc1,'Value',0)
    end
    if FEATURE_VECTOR.pc2
        set(handles.checkbox_pc2,'Value',1)
    else
        set(handles.checkbox_pc2,'Value',0)
    end
    if FEATURE_VECTOR.pc3
        set(handles.checkbox_pc3,'Value',1)
    else
        set(handles.checkbox_pc3,'Value',0)
    end
    if FEATURE_VECTOR.pc4
        set(handles.checkbox_pc4,'Value',1)
    else
        set(handles.checkbox_pc4,'Value',0)
    end
    if FEATURE_VECTOR.pce1
        set(handles.checkbox_pce1,'Value',1)
    else
        set(handles.checkbox_pce1,'Value',0)
    end
    if FEATURE_VECTOR.pce2
        set(handles.checkbox_pce2,'Value',1)
    else
        set(handles.checkbox_pce2,'Value',0)
    end
    if FEATURE_VECTOR.pce3
        set(handles.checkbox_pce3,'Value',1)
    else
        set(handles.checkbox_pce3,'Value',0)
    end
    if FEATURE_VECTOR.pce4
        set(handles.checkbox_pce4,'Value',1)
    else
        set(handles.checkbox_pce4,'Value',0)
    end
    if FEATURE_VECTOR.pca1
        set(handles.checkbox_pca1,'Value',1)
    else
        set(handles.checkbox_pca1,'Value',0)
    end
    if FEATURE_VECTOR.pca2
        set(handles.checkbox_pca2,'Value',1)
    else
        set(handles.checkbox_pca2,'Value',0)
    end
    if FEATURE_VECTOR.pca3
        set(handles.checkbox_pca3,'Value',1)
    else
        set(handles.checkbox_pca3,'Value',0)
    end
    if FEATURE_VECTOR.pca4
        set(handles.checkbox_pca4,'Value',1)
    else
        set(handles.checkbox_pca4,'Value',0)
    end
    if FEATURE_VECTOR.peak
        set(handles.checkbox_peak,'Value',1)
    else
        set(handles.checkbox_peak,'Value',0)
    end
    if FEATURE_VECTOR.valley
        set(handles.checkbox_valley,'Value',1)
    else
        set(handles.checkbox_valley,'Value',0)
    end
    if FEATURE_VECTOR.amp
        set(handles.checkbox_amp,'Value',1)
    else
        set(handles.checkbox_amp,'Value',0)
    end
    if FEATURE_VECTOR.energy
        set(handles.checkbox_e,'Value',1)
    else
        set(handles.checkbox_e,'Value',0)
    end
    if FEATURE_VECTOR.area
        set(handles.checkbox_a,'Value',1)
    else
        set(handles.checkbox_a,'Value',0)
    end
    if FEATURE_VECTOR.params
        set(handles.checkbox_params,'Value',1)
    else
        set(handles.checkbox_params,'Value',0)
    end
    if FEATURE_VECTOR.slice
        set(handles.checkbox_slice,'Value',1)
    else
        set(handles.checkbox_slice,'Value',0)
    end
end

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------
function varargout = czfeaturespace_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------
function checkbox_pc1_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pc2_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pc3_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pc4_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pce1_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pce2_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pce3_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pce4_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pca1_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pca2_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pca3_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_pca4_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_peak_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_e_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_valley_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_amp_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_a_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_params_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function checkbox_slice_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function pushbutton_ok_Callback(hObject, eventdata, handles)

% Create output
feature_names = {};
feature_vector = struct();
feature_vector.pc1 = get(handles.checkbox_pc1,'Value');
if feature_vector.pc1
    feature_names = [feature_names {'EL1 PC1'} {'EL2 PC1'} {'EL3 PC1'} {'EL4 PC1'}];
end
feature_vector.pc2 = get(handles.checkbox_pc2,'Value');
if feature_vector.pc2
    feature_names = [feature_names {'EL1 PC2'} {'EL2 PC2'} {'EL3 PC2'} {'EL4 PC2'}];
end
feature_vector.pc3 = get(handles.checkbox_pc3,'Value');
if feature_vector.pc3
    feature_names = [feature_names {'EL1 PC3'} {'EL2 PC3'} {'EL3 PC3'} {'EL4 PC3'}];
end
feature_vector.pc4 = get(handles.checkbox_pc4,'Value');
if feature_vector.pc4
    feature_names = [feature_names {'EL1 PC4'} {'EL2 PC4'} {'EL3 PC4'} {'EL4 PC4'}];
end

feature_vector.pce1 = get(handles.checkbox_pce1,'Value');
if feature_vector.pce1
    feature_names = [feature_names {'EL1 PCE1'} {'EL2 PCE1'} {'EL3 PCE1'} {'EL4 PCE1'}];
end
feature_vector.pce2 = get(handles.checkbox_pce2,'Value');
if feature_vector.pce2
    feature_names = [feature_names {'EL1 PCE2'} {'EL2 PCE2'} {'EL3 PCE2'} {'EL4 PCE2'}];
end
feature_vector.pce3 = get(handles.checkbox_pce3,'Value');
if feature_vector.pce3
    feature_names = [feature_names {'EL1 PCE3'} {'EL2 PCE3'} {'EL3 PCE3'} {'EL4 PCE3'}];
end
feature_vector.pce4 = get(handles.checkbox_pce4,'Value');
if feature_vector.pce4
    feature_names = [feature_names {'EL1 PCE4'} {'EL2 PCE4'} {'EL3 PCE4'} {'EL4 PCE4'}];
end

feature_vector.pca1 = get(handles.checkbox_pca1,'Value');
if feature_vector.pca1
    feature_names = [feature_names {'EL1 PCA1'} {'EL2 PCA1'} {'EL3 PCA1'} {'EL4 PCA1'}];
end
feature_vector.pca2 = get(handles.checkbox_pca2,'Value');
if feature_vector.pca2
    feature_names = [feature_names {'EL1 PCA2'} {'EL2 PCA2'} {'EL3 PCA2'} {'EL4 PCA2'}];
end
feature_vector.pca3 = get(handles.checkbox_pca3,'Value');
if feature_vector.pca3
    feature_names = [feature_names {'EL1 PCA3'} {'EL2 PCA3'} {'EL3 PCA3'} {'EL4 PCA3'}];
end
feature_vector.pca4 = get(handles.checkbox_pca4,'Value');
if feature_vector.pca4
    feature_names = [feature_names {'EL1 PCA4'} {'EL2 PCA4'} {'EL3 PCA4'} {'EL4 PCA4'}];
end

feature_vector.peak = get(handles.checkbox_peak,'Value');
if feature_vector.peak
    feature_names = [feature_names {'EL1 PEAK'} {'EL2 PEAK'} {'EL3 PEAK'} {'EL4 PEAK'}];
end
feature_vector.valley = get(handles.checkbox_valley,'Value');
if feature_vector.valley
    feature_names = [feature_names {'EL1 VALLEY'} {'EL2 VALLEY'} {'EL3 VALLEY'} {'EL4 VALLEY'}];
end
feature_vector.amp = get(handles.checkbox_amp,'Value');
if feature_vector.amp
    feature_names = [feature_names {'EL1 AMP'} {'EL2 AMP'} {'EL3 AMP'} {'EL4 AMP'}];
end
feature_vector.slice = get(handles.checkbox_slice,'Value');
if feature_vector.slice
    feature_names = [feature_names {'EL1 SLICE'} {'EL2 SLICE'} {'EL3 SLICE'} {'EL4 SLICE'}];
end

feature_vector.energy = get(handles.checkbox_e,'Value');
if feature_vector.energy
    feature_names = [feature_names {'EL1 E'} {'EL2 E'} {'EL3 E'} {'EL4 E'}];
end
feature_vector.area = get(handles.checkbox_a,'Value');
if feature_vector.area
    feature_names = [feature_names {'EL1 A'} {'EL2 A'} {'EL3 A'} {'EL4 A'}];
end
feature_vector.params = get(handles.checkbox_params,'Value');
if feature_vector.params
    feature_names = [feature_names {'EL1 P1'} {'EL2 P1'} {'EL3 P1'} {'EL4 P1'} ...
        {'EL1 P2'} {'EL2 P2'} {'EL3 P2'} {'EL4 P2'}];
end

% Assign output
global FEATURE_NAMES
FEATURE_NAMES = feature_names;
global FEATURE_VECTOR
FEATURE_VECTOR = feature_vector;

% Close GUI
delete(handles.figure1)

% -------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)

% Close GUI
delete(handles.figure1)