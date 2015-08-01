function varargout = agannot_gui(varargin)
% AGANNOT_GUI M-file for agannot_gui.fig
%      AGANNOT_GUI, by itself, creates a new AGANNOT_GUI or raises the existing
%      singleton*.
%
%      H = AGANNOT_GUI returns the handle to a new AGANNOT_GUI or the handle to
%      the existing singleton*.
%
%      AGANNOT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGANNOT_GUI.M with the given input arguments.
%
%      AGANNOT_GUI('Property','Value',...) creates a new AGANNOT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before agannot_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to agannot_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%AGANNOT_GUI   Data annotation for ANALYSE_SNIFFWHISK GUI.
%   AGANNOT_GUI enables the creation of annotation objects for
%   ANALYSE_SNIFFWHISK GUI. The objects contain a start point, an end point
%   and a string. These parameters can be  set interactively. Start point
%   and end point can be set by pressing 's' and 'e' keys in the parrent
%   GUI. The annotation object is passed back to the parrent GUI.
%
%   See also ANALYSE_SNIFFWHISK and AGANNOT.

% Edit the above text to modify the response to help agannot_gui

% Last Modified by GUIDE v2.5 16-Feb-2011 10:43:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @agannot_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @agannot_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    if ~isempty(varargin)   % set GUI position
        pd = varargin{1};
        ps = get(varargout{1},'Position');
        lf = pd(1) + pd(3) / 2 - ps(3) / 2;
        bt = pd(2) + pd(4) / 2 - ps(4) / 2;
        set(varargout{1},'Position',[lf bt ps(3) ps(4)]);
    end
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% -------------------------------------------------------------------------
function agannot_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% Executes just before agannot_gui is made visible.
% This function has no output args, see OutputFcn.

% Choose default command line output for agannot_gui
handles.output = hObject;

% Parrent GUI
if ~isempty(varargin)
    handles.parrentgui = varargin{2};
else
    handles.parrentgui = 0;
end

% Previous settings
if length(varargin) > 2
    set(handles.edit1,'String',num2str(varargin{3}))
    set(handles.edit2,'String',num2str(varargin{4}))
    set(handles.edit3,'String',num2str(varargin{5}))
end

% Update handles structure
guidata(hObject, handles);



% -------------------------------------------------------------------------
function varargout = agannot_gui_OutputFcn(hObject, eventdata, handles)
% Outputs from this function are returned to the command line. 

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pushbutton Callbacks                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Start Point
% -------------------------------------------------------------------------
function pushbutton1_Callback(hObject, eventdata, handles)

% Wait for 's' keypress
w = 1;
on = get(handles.parrentgui,'Name');
set(handles.parrentgui,'Name','Press ''s'' to determine start point')
while w
    uiwait(handles.parrentgui)
    ch = get(handles.parrentgui,'CurrentCharacter');
    if ch == 's'
        w = 0;
    end
end
set(handles.parrentgui,'Name',on)

% Set start point
cp = get(get(handles.parrentgui,'CurrentAxes'),'CurrentPoint');
set(handles.edit1,'String',num2str(cp(1)))

% Return control to annotation GUI
figure(handles.figure1)



% -------------------------------------------------------------------------
% End Point
% -------------------------------------------------------------------------
function pushbutton2_Callback(hObject, eventdata, handles)

% Wait for 'e' keypress
w = 1;
on = get(handles.parrentgui,'Name');
set(handles.parrentgui,'Name','Press ''e'' to determine end point')
while w
    uiwait(handles.parrentgui)
    ch = get(handles.parrentgui,'CurrentCharacter');
    if ch == 'e'
        w = 0;
    end
end
set(handles.parrentgui,'Name',on)

% Set end point
cp = get(get(handles.parrentgui,'CurrentAxes'),'CurrentPoint');
set(handles.edit2,'String',num2str(cp(1)))

% Return control to annotation GUI
figure(handles.figure1)



% -------------------------------------------------------------------------
% OK
% -------------------------------------------------------------------------
function pushbutton3_Callback(hObject, eventdata, handles)

% Output
psp = get(handles.edit1,'String');
if ~isempty(psp)
    sp = str2double(psp);
else
    sp = NaN;
end
pep = get(handles.edit2,'String');
if ~isempty(pep)
    ep = str2double(pep);
else
    ep = NaN;
end
str = get(handles.edit3,'String');

% Delete GUI
figure1_CloseRequestFcn(hObject, eventdata, handles)

% Assign output
global ANNOTOUT
annotout.sp = sp;
annotout.ep = ep;
annotout.str = str;
ANNOTOUT = annotout;



% -------------------------------------------------------------------------
% CANCEL
% -------------------------------------------------------------------------
function pushbutton4_Callback(hObject, eventdata, handles)

% Close GUI
figure1_CloseRequestFcn(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Edit Text Callbacks                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit1_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit1_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function edit2_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit2_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function edit3_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
function edit3_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GUI Figure Callbacks                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Close Request Function
% -------------------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)

% Clear GUI
set(handles.edit1,'String','')
set(handles.edit2,'String','')
set(handles.edit3,'String','')

% Assign output
global ANNOTOUT
annotout.sp = NaN;
annotout.ep = NaN;
annotout.str = '';
ANNOTOUT = annotout;

% Delete GUI
figure1_DeleteFcn(hObject, eventdata, handles)



% -------------------------------------------------------------------------
% Delete Function
% -------------------------------------------------------------------------
function figure1_DeleteFcn(hObject, eventdata, handles)

% Delete GUI
delete(handles.figure1)