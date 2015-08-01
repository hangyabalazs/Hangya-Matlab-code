function varargout = b_ica_layouteditor(varargin)
% B_ICA_LAYOUTEDITOR Application M-file for b_ica_layouteditor.fig
%    FIG = B_ICA_LAYOUTEDITOR launch b_ica_layouteditor GUI.
%    B_ICA_LAYOUTEDITOR('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 31-Aug-2004 15:05:25

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    
    % Load last setting
    global DATAPATH
    try
        ff = fullfile(DATAPATH,'ICA\ica_gui4\settings\gui_layout.mat');
    end
    if b_isfilename(ff)
        load(ff)
        nrows = gui_layout{1};
        ncols = gui_layout{2};
        lenr = length(nrows);
        lenc = length(ncols);
        set(handles.nrows,'String',nrows)
        set(handles.ncols,'String',ncols)
        if strcmp(layout_type,'rowtype')
            set(handles.rowtype,'Value',1);
            set(handles.columntype,'Value',0);
            if lenc > 1
                set(handles.rowtype,'Enable','off');
                set(handles.columntype,'Enable','off');
            elseif lenr == 1 & lenc == 1
                set(handles.rowtype,'Enable','on');
                set(handles.columntype,'Enable','on');
            end
        elseif strcmp(layout_type,'columntype')
            set(handles.rowtype,'Value',0);
            set(handles.columntype,'Value',1);
            if lenr > 1
                set(handles.rowtype,'Enable','off');
                set(handles.columntype,'Enable','off');
            elseif lenr == 1 & lenc == 1
                set(handles.rowtype,'Enable','on');
                set(handles.columntype,'Enable','on');
            end
        end
    else
        warndlg('Cannot load layout settings: permission denied or file does not exist. Setting values to default.',...
            'Loading failure');
        default_Callback(fig,[],handles,varargin)
    end

	if nargout > 0
		varargout{1} = fig;
	end

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
% 'NUMBER OF ROWS' edit text callback
% --------------------------------------------------------------------
function varargout = nrows_Callback(h,eventdata,handles,varargin)

% Set radiobuttons' properties
nrows = str2num(get(handles.nrows,'String'));
ncols = str2num(get(handles.ncols,'String'));
lenr = length(nrows);
lenc = length(ncols);
set_radiobuttons(h,handles,lenr,lenc)



% --------------------------------------------------------------------
% 'NUMBER OF COLUMNS' edit text callback
% --------------------------------------------------------------------
function varargout = ncols_Callback(h,eventdata,handles,varargin)

% Set radiobuttons' properties
nrows = str2num(get(handles.nrows,'String'));
ncols = str2num(get(handles.ncols,'String'));
lenr = length(nrows);
lenc = length(ncols);
set_radiobuttons(h,handles,lenr,lenc)



% --------------------------------------------------------------------
% 'ROWTYPE' radiobutton callback
% --------------------------------------------------------------------
function varargout = rowtype_Callback(h,eventdata,handles,varargin)
set(handles.columntype,'Value',0);



% --------------------------------------------------------------------
% 'COLUMNTYPE' radiobutton callback
% --------------------------------------------------------------------
function varargout = columntype_Callback(h,eventdata,handles,varargin)
set(handles.rowtype,'Value',0);




% --------------------------------------------------------------------
% DEFAULT button callback - reset values
% --------------------------------------------------------------------
function varargout = default_Callback(h,eventdata,handles,varargin)

% Reset default values
set(handles.nrows,'String','6');
set(handles.ncols,'String','1 1 1 2 2 2')
set(handles.rowtype,'Value',1)
set(handles.columntype,'Value',0)

% Set radiobuttons' properties
set(handles.rowtype,'Value',1);
set(handles.rowtype,'Enable','off');
set(handles.columntype,'Value',0);
set(handles.columntype,'Enable','off');



% --------------------------------------------------------------------
% OK button callback - save setings and delete GUI
% --------------------------------------------------------------------
function varargout = ok_Callback(h,eventdata,handles,varargin)

% Get edit text strings and check if they have been set properly
nrows = str2num(get(handles.nrows,'String'));
ncols = str2num(get(handles.ncols,'String'));
lenr = length(nrows);
lenc = length(ncols);
if (lenc > 1 & lenr ~= 1) | (lenr > 1 & lenc ~= 1)
    warndlg('Invalid settings: at least one of the inputs has to be scalar.','Warning');
    return
end
if (lenr == 1 & ~isequal(lenc,1) & ~isequal(nrows,lenc)) | (lenc == 1 &~isequal(lenr,1) & ~isequal(ncols,lenr))
    warndlg('Invalid  settings: dimension mismatch.','Warning')
    return
end
isrt = get(handles.rowtype,'Value');
if isrt
    layout_type = 'rowtype';
else
    layout_type = 'columntype';
end

% Change layout
global ICA_GUI_HANDLE
h2 = ICA_GUI_HANDLE;
figure(h2)       % set current figure
handles2 = guidata(h2);
A = handles2.axes;
b_ica_gui4('layout',h2,handles2,nrows,ncols,handles2.axes)
handles2 = guidata(h2);       % set current axes
A = handles2.axes;
axes(A(1))
set(A(1),'XColor','red')     % set box color
set(A(1),'YColor','red')

% Move legends
for x = 1:length(A)
    [L,legendpos] = b_ica_gui4('legendcheck',handles2.figure1,A(x));
    if ~isempty(L)
        old_axes_units = get(A(x),'Units');
        old_legend_units = get(L,'Units');
        set(A(x),'Units','points')
        set(L,'Units','points')
        pos = get(L,'Position');
        width = pos(3);
        hight = pos(4);
        lpos = b_ica_gui4('legendposition',A(x),legendpos,width,hight);
        set(L,'Position',lpos,'Units',old_legend_units)
        set(A(x),'Units',old_axes_units)
    end
end

% Move buttons
button_handles = findobj(get(handles2.figure1,'Children'),'Style','pushbutton');
button_handles = flipud(button_handles);
if ~isempty(button_handles)     % remember old 'Units'
    old_units_buttons = get(button_handles(1),'Units');
end
old_axes_units = get(A(x),'Units');
set(button_handles,'Units','pixels');   % set 'Units'
set(A,'Units','pixels');

next = 1;
ax = 0;
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    if ~isequal(usd.axes_handle,ax)
        next = 1;
    end
    ax = usd.axes_handle;
    lenb = usd.button_number;
    ax_pos = get(ax,'Position');
    width = ax_pos(3) / lenb;
    left = ax_pos(1) + (next - 1) * width;
    buttom = ax_pos(2) + ax_pos(4);
    hight = 20;
    buttons_pos = [left buttom width hight];
    set(button_handles(i),'Position',buttons_pos);
    next = next + 1;
end

if ~isempty(button_handles)
    set(button_handles,'Units',old_units_buttons);   % reset 'Units'
end
set(A,'Units',old_axes_units)

% Save settings and delete GUI
nrows = get(handles.nrows,'String');
ncols = get(handles.ncols,'String');
gui_layout = [{nrows} {ncols}];
global DATAPATH
ff = fullfile(DATAPATH,'ICA\ica_gui4\settings\gui_layout.mat');
str = ['save ' ff ' gui_layout layout_type'];
eval(str)
delete(handles.figure1)



% --------------------------------------------------------------------
% CANCEL button callback - delete GUI
% --------------------------------------------------------------------
function varargout = cancel_Callback(h,eventdata,handles,varargin)

% Delete GUI
delete(handles.figure1)



% --------------------------------------------------------------------
% Set radiobutton properties
% --------------------------------------------------------------------
function set_radiobuttons(h,handles,lenr,lenc)

% Set radiobuttons' 'Value' and 'Enable' properties
if lenr > 1 & lenc > 1
    set(handles.rowtype,'Value',0);
    set(handles.columntype,'Value',0);
    set(handles.rowtype,'Enable','off');
    set(handles.columntype,'Enable','off');
elseif lenr > 1 & lenc == 1
    set(handles.rowtype,'Value',0);
    set(handles.columntype,'Value',1);
    set(handles.rowtype,'Enable','off');
    set(handles.columntype,'Enable','off');
elseif lenr == 1 & lenc > 1
    set(handles.rowtype,'Value',1);
    set(handles.columntype,'Value',0);
    set(handles.rowtype,'Enable','off');
    set(handles.columntype,'Enable','off');
elseif lenr == 1 & lenc == 1
    set(handles.rowtype,'Enable','on');
    set(handles.columntype,'Enable','on');
end