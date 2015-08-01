function varargout = b_empty_gui2004(varargin)
%EMPTY_GUI2004   Batch processing on hippocampus data files.
%   You are able to specify the following parameters via EMPTY_GUI2004:
%       Running function: enter the name of the running function. See
%                         EMPTY_RUNFILE2004 for an example.
%       Input arguments:  white space sepatated list. Elements should be 
%                         names of base workspace variables.
%       Text file name:   name of text output. Leave it empty if no text
%                         output is required.
%       Directory containing the data: data directory.
%       Directory for the results: results' directory.
%       Interval set:     auto - get interval limits from datinx files
%                         manually - inport them from Command Window.
%       Name of files to save: saved file names are generated automatically
%                         from this input and the file identifiers
%
%   EMPTY_GUI2004 calls EMPTY2004.
%
%   See also EMPTY2004, EMPTY_RUNFILE2004 and EMPTY_NEW.


% Last Modified by GUIDE v2.0 19-Apr-2004 16:46:25

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
    
    % Set IntervalSet to 'automatic'
    set(handles.radiobutton2,'Value',1)
    
    % Preallocate some fields
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
fstr = findstr(handles.InputArg,',');
fstr(end+1) = length(handles.InputArg) + 1;
lfstr = length(fstr);
ias = cell(1,lfstr);
ia = cell(1,lfstr);
next = 1;
for i = 1:lfstr
    ias{i} = handles.InputArg(next:fstr(i)-1);
    next = fstr(i) + 1;
    cmnd = ['assignin(''caller'',''pia'',' ias{i} ');'];
    evalin('base',cmnd)
    ia{i} = pia;
end
handles.InputArg = ia;
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for editing DIRECTORY containing the DATA
% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
handles.DataDir = get(handles.edit2,'String');
if ~(handles.DataDir(1)==''''&handles.DataDir(end)=='''')
    cmnd = ['assignin(''caller'',''dd'',' handles.DataDir ');'];
    evalin('base',cmnd)
    handles.DataDir = dd;
end
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for editing DIRECTORY for the RESULTS
% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
handles.ResultDir = get(handles.edit3,'String');
if ~(handles.ResultDir(1)==''''&handles.ResultDir(end)=='''')
    cmnd = ['assignin(''caller'',''rd'',' handles.ResultDir ');'];
    evalin('base',cmnd)
    handles.ResultDir = rd;
end
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
% Callback for TEXT OUTPUT set 
% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
handles.TextOutput = get(handles.edit4,'String');
if ~(handles.TextOutput(1)==''''&handles.TextOutput(end)=='''')
    cmnd = ['assignin(''caller'',''to'',' handles.TextOutput ');'];
    evalin('base',cmnd)
    handles.TextOutput = to;
end
guidata(h,handles);

% --------------------------------------------------------------------
% Callback for START button - run
% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

int = get(handles.radiobutton1,'Value');
if int
    IntervalSet = 'manual';
else
    IntervalSet = 'auto';
end

b_empty2004(handles.RunFile,handles.DataDir,handles.ResultDir,IntervalSet,...
        handles.TextOutput,handles.Name,handles.InputArg);

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


% --------------------------------------------------------------------
% Callback for editing NAME OF FILES TO SAVE
% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)
handles.Name = get(handles.edit6,'String');
guidata(h,handles)