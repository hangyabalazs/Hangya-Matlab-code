function varargout = b_theta_gui(varargin)
%THETA_GUI  Display results of theta selection.
%   THETA_GUI plots theta-delta ratio, distribution of theta-delta ratio and wavelet
%   with maximum locations (output of THETASELECTOR3 - see THETASELECTOR3 for details).
%   The functions loads the results of THETASELECTORRUN.
%
%   For more details on theta selection, see help for the related functions (below).
%
%   See also THETASELECTOR3, THETASELECTOR_BETA3, THETASELECTORRUN and THETA.

% Last Modified by GUIDE v2.0 07-Apr-2004 13:35:53

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
    
    eventdata = [];
    setliststring(fig, eventdata, handles, varargin);
    
    % Store default GUI position in the handles structure for the Resize Function
    old_units = get(handles.figure1,'Units');
    set(handles.figure1,'Units','pixels');
    default_gui_pos = get(handles.figure1,'Position');
    set(handles.figure1,'Units',old_units);
    handles.gui_pos = default_gui_pos;    
    guidata(fig,handles);
    
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
        errordlg(lasterr,'Error message');
	end

end




% --------------------------------------------------------------------
% LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)
if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
    lb(h, eventdata, handles, varargin)
end

% --------------------------------------------------------------------
function varargout = lb(h, eventdata, handles, varargin)

% Load OM
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');

fln = file_list{index_selected};

ff = fullfile(DATAPATH,'Wavelet\thetaselection\matrix',file_list{index_selected});
load(ff)

% Plot on axis1 - theta power / delta power
axes(handles.axes1)
bdf1 = get(handles.axes1,'ButtonDownFcn');

ratio = Out(3,:);
plot(ratio,'r')
y_lim = ylim;
x_lim = xlim;
axis([x_lim(1) length(ratio) y_lim(1) y_lim(2)]);
tt = '\theta power / \delta power';
xcoord = 1 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 2;
ycoord = 3 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
text(xcoord,ycoord,tt,'Color',[0 0 0]);

% Plot on axis2 - distribution of theta power / delta power values
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');

[x y] = hist(ratio,100);
bar(y,x);
y_lim = ylim;
x_lim = xlim;
tt = 'Distribution of \theta power / \delta power';
xcoord = 1 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 2;
ycoord = 3 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
text(xcoord,ycoord,tt,'Color',[0 0 0]);

% Plot on axis3 - wavelet
axes(handles.axes3)
fln2 = [fln(1:end-4) '.jpg'];
ff = fullfile(DATAPATH,'Wavelet\thetaselection\jpg',fln2);
str = ['A = imread(''' ff ''');'];
eval(str)
imagesc(A)
old_units_axes = get(handles.axes3,'Units');
set(handles.axes3,'Units','normalized');
axis([156 1088 68 821])
set(handles.axes3,'Units',old_units_axes);
y_lim = ylim;
x_lim = xlim;
tt = 'EEG Wavelet';
xcoord = 1 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 2;
ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
text(xcoord,ycoord,tt,'Color',[1 1 1]);

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
set(handles.axes1,'ButtonDownFcn',bdf1);

% Set keypress function
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

guidata(h,handles);

% --------------------------------------------------------------------
% RESIZE function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
gui_pos = handles.gui_pos;

% Remember the old Units properties
old_units_fig = get(handles.figure1,'Units');
old_units_text = get(handles.text1,'Units');
old_units_axes = get(handles.axes1,'Units');
old_units_listbox = get(handles.listbox1,'Units');

% Set Units to pixels
set(handles.figure1,'Units','pixels');
pos = get(handles.figure1,'Position');

% Compensate if too narrow
if pos(3) < 300
    if pos(1) == gui_pos(1)
        pos(3) = 300;
        set(handles.figure1,'Position',pos);
    else
        pos(1) = pos(1) + pos(3) - 300;
        pos(3) = 300;
        set(handles.figure1,'Position',pos);
    end
end

% Compensate if too low
if pos(4) < 455
    if pos(2) == gui_pos(2)
        pos(4) = 455;
        set(handles.figure1,'Position',pos);
    else
        pos(2) = pos(2) + pos (4) - 455;
        pos(4) = 455;
        set(handles.figure1,'Position',pos);
    end
end

% New position of the listbox
pos_list = cell(1,1);
for t = 1:1
    eval(['set(handles.listbox',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_list{t} = get(handles.listbox',int2str(t),',''Position'');']);
    eval(['set(handles.listbox',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_list{t}(1) pos(4)-gui_pos(4)+pos_list{t}(2) pos_list{t}(3) pos_list{t}(4)]);']);
end

% New position of the static text above the listbox
pos_text = cell(1,3);
set(handles.text1,'Units','pixels');
pos_text{1} = get(handles.text1,'Position');
set(handles.text1,...
    'Position',[pos(3)-gui_pos(3)+pos_text{1}(1) pos(4)-gui_pos(4)+pos_text{1}(2) pos_text{1}(3) pos_text{1}(4)]);

% New position of the axes
pos_axes = cell(1,3);
set(handles.axes1,'Units','pixels');
pos_axes{1} = get(handles.axes1,'Position');
set(handles.axes1,...
    'Position',[pos_axes{1}(1) pos_axes{1}(2)+((pos(4)-gui_pos(4))/2)...
        pos(3)-gui_pos(3)+pos_axes{1}(3) ((pos(4)-gui_pos(4))/2)+pos_axes{1}(4)]);
set(handles.axes2,'Units','pixels');
pos_axes{2} = get(handles.axes2,'Position');
set(handles.axes2,...
    'Position',[pos_axes{2}(1) pos_axes{2}(2)...
        pos(3)-gui_pos(3)+pos_axes{2}(3) ((pos(4)-gui_pos(4))/2)+pos_axes{2}(4)]);
set(handles.axes3,'Units','pixels');
pos_axes{3} = get(handles.axes3,'Position');
set(handles.axes3,...
    'Position',[pos_axes{3}(1) pos_axes{3}(2)...
        pos(3)-gui_pos(3)+pos_axes{3}(3) pos_axes{3}(4)]);

% New position of the static texts above the axes
% set(handles.text2,'Units','pixels');
% pos_text{2} = get(handles.text2,'Position');
% set(handles.text2,...
%     'Position',[((pos(3)-gui_pos(3))/2)+pos_text{2}(1) pos(4)-gui_pos(4)+pos_text{2}(2)...
%         pos_text{2}(3) pos_text{2}(4)]);
% set(handles.text1,'Units','pixels');
% pos_text{1} = get(handles.text1,'Position');
% set(handles.text1,...
%     'Position',[((pos(3)-gui_pos(3))/2)+pos_text{1}(1) ((pos(4)-gui_pos(4))/2)+pos_text{1}(2)...
%         pos_text{1}(3) pos_text{1}(4)]);

% Reposition GUI on screen
% movegui(h,'onscreen')

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);
set(handles.text1,'Units',old_units_text);
set(handles.axes1,'Units',old_units_axes);
set(handles.listbox1,'Units',old_units_listbox);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);

% --------------------------------------------------------------------
% AXES1 BUTTONDOWN function - step one cell
% --------------------------------------------------------------------
function varargout = axes1_ButtonDownFcn(h, eventdata, handles, varargin)
index_selected = get(handles.listbox1,'Value');
list = get(handles.listbox1,'String');
lenlist = length(list);
seltyp = get(handles.figure1,'SelectionType');
switch seltyp
case {'normal','open'}
    if index_selected + 1 <= lenlist
        set(handles.listbox1,'Value',index_selected+1)
        lb(h, eventdata, handles, varargin);
    end
otherwise
    if index_selected - 1 >= 1
        set(handles.listbox1,'Value',index_selected-1)
        lb(h, eventdata, handles, varargin);
    end
end

% --------------------------------------------------------------------
% AXES2 BUTTONDOWN function - step one figure
% --------------------------------------------------------------------
function varargout = axes2_ButtonDownFcn(h, eventdata, handles, varargin)
index_selected = get(handles.listbox2,'Value');
list = get(handles.listbox2,'String');
lenlist = length(list);
seltyp = get(handles.figure1,'SelectionType');
switch seltyp
case {'normal','open'}
    if index_selected + 1 <= lenlist
        set(handles.listbox2,'Value',index_selected+1)
        lb2(h, eventdata, handles, varargin);
    end
otherwise
    if index_selected - 1 >= 1
        set(handles.listbox2,'Value',index_selected-1)
        lb2(h, eventdata, handles, varargin);
    end
end
    
% ---------------------------------------------------------------------
% SETLISTSTRING subfunction - set string of listbox1
% ---------------------------------------------------------------------
function varargout = setliststring(h, eventdata, handles, varargin)
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
ff = fullfile(DATAPATH,'Wavelet\thetaselection\matrix');
k = dir(ff);
lenlist = length(k)-2;
liststring = {};
for t = 1:lenlist
    if ~k(t+2).isdir
        liststring{end+1} = k(t+2).name;
    end
end
set(handles.listbox1,'String',liststring);
guidata(h,handles)