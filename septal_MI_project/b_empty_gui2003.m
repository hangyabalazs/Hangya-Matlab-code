function varargout = b_empty_gui2003(varargin)
%EMPTY_GUI2003   Batch processing on hippocampus data files.
%   You are able to specify the following parameters via EMPTY_GUI2003:
%       Running function: enter the name of the running function. You
%                         can create a running function from EMPTY2003.m.
%       Input arguments:  white space sepatated list. Elements can be
%                         either numeric values or names of base workspace
%                         variables.
%       Text file name:   name of text output. Leave it empty if no text
%                         output is required.
%       Directory containing the data: data directory.
%       Directory for the results: results' directory.
%       Interval set:     auto - get interval limits from datinx files
%                         manually - inport them from Command Window.
%       Evoked theta, spontanous theta, no theta: check if examination
%                         of the specified segment is required
%
%   See also EMPTY2003, EMPTY_GUI2004 and EMPTY_NEW.

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
    
    handles.TextOutput = [];
    handles.InputArg = [];
    guidata(fig, handles);

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


% --------------------------------------------------------------------
% Callback for editing INPUT ARGUMENTS
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
handles.InputArg = get(handles.edit1,'String');
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for editing DIRECTORY containing the DATA
% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
handles.DataDir = get(handles.edit2,'String');
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for editing DIRECTORY for the RESULTS
% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
handles.ResultDir = get(handles.edit3,'String');
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for MANUAL INTERVAL SET 
% --------------------------------------------------------------------
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)
set(handles.radiobutton2,'Value',0)

% --------------------------------------------------------------------
% Callback for AUTOMATIC INTERVAL SET 
% --------------------------------------------------------------------
function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)
set(handles.radiobutton1,'Value',0)

% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
% Callback for TEXT OUTPUT set 
% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
handles.TextOutput = get(handles.edit4,'String');
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for START button - run
% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
EvTheta = get(handles.checkbox1,'Value');
SpTheta = get(handles.checkbox2,'Value');
NoTheta = get(handles.checkbox3,'Value');
EvSpNo = [EvTheta SpTheta NoTheta];

int = get(handles.radiobutton1,'Value');
if int
    IntervalSet = 'manual';
else
    IntervalSet = 'auto';
end

str = [handles.RunFile '(handles.DataDir,handles.ResultDir,IntervalSet,'...
        'handles.TextOutput,EvSpNo,handles.InputArg)'];
eval(str);


% --------------------------------------------------------------------
% Callback for CANCEL button - deletes GUI
% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
delete(handles.figure1)

% --------------------------------------------------------------------
% Callback for editing name of RUNNING FUNCTION
% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)
handles.RunFile = get(handles.edit5,'String');
guidata(h,handles)