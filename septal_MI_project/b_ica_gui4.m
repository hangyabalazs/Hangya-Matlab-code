function varargout = b_ica_gui4(varargin)
% B_ICA_GUI4 Application M-file for b_ica_gui4.fig
%    FIG = B_ICA_GUI4 launch b_ica_gui4 GUI.
%    B_ICA_GUI4('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 08-Sep-2004 17:53:29



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       L A U N C H   G U I                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	% Export GUI handle ('fig')
    global ICA_GUI_HANDLE
    ICA_GUI_HANDLE = fig;
    
    % Store default GUI position in the handles structure for the Resize Function
    old_units = get(handles.figure1,'Units');
    set(handles.figure1,'Units','pixels');
    default_gui_pos = get(handles.figure1,'Position');
    set(handles.figure1,'Units',old_units);
    handles.gui_pos = default_gui_pos;    
    guidata(fig,handles);

    % LAYOUT
    nrows = [6];
    ncols = [1 1 1 2 2 2];
    handles.type = 'rowtype';
    guidata(fig,handles)
    layout(fig,handles,nrows,ncols)
    handles = guidata(fig);
    
    % Set current axes
    A = handles.axes;
    axes(A(1))
    set(A(1),'XColor','red')     % set box color
    set(A(1),'YColor','red')
    
    % UIMENUS
    menus(fig,handles)
    
    % Refresh handles structure 
    handles = guidata(fig);
    
    % Set default clustering method to 'ward'
    handles.clustering_method = 'ward';
    
    % Initialize 'figure' field
    handles.figure = [];
    
    % Initialize 'display' fields
    handles.guidisplay = struct('name',[],'settings',{},'button',[],'stringinput',[],'keypress',...
        [],'opening_button',[]);
    emp.name = '';
    emp.settings = {};
    emp.button = [];
    emp.stringinput = '';
    emp.keypress = '';
    emp.opening_button = [];
    if length(ncols) > 1
        axno = sum(ncols);
    elseif length(nrows) > 1
        axno = sum(nrows);
    else
        axno = ncols * nrows;
    end
    for i = 1:axno
        handles.guidisplay(i) = emp;
    end
    handles.figdisplay = struct('name',[],'settings',{},'button',[],'stringinput',[],'keypress',...
        [],'opening_button',[]);

    % Initialize 'greyplot' fields, 'line_handles' field
    handles.greyplot1_handle = [];
    handles.greyplot2_handle = [];
    handles.greyplot3_handle = [];
    handles.greyplot4_handle = [];
    handles.line_handles = [];
    guidata(fig,handles)
    
    
    
    
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
        errordlg(lasterr,'Error message');
	end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           L A Y O U T                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------------------------------------------------------------
% RESIZE function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h,eventdata,handles,varargin)

% Get default GUI position
gui_pos = handles.gui_pos;

% Get pushbutton handles
button_handles = findobj(get(handles.figure1,'Children'),'Style','pushbutton');
button_handles = flipud(button_handles);

% Remember the old 'Units' properties
old_units_fig = get(handles.figure1,'Units');
if ~isempty(button_handles)
    old_units_buttons = get(button_handles(1),'Units');
end

% Set 'Units' to pixels
set(handles.figure1,'Units','pixels');
set(button_handles,'Units','pixels');
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
if pos(4) < 300
    if pos(2) == gui_pos(2)
        pos(4) = 300;
        set(handles.figure1,'Position',pos);
    else
        pos(2) = pos(2) + pos (4) - 300;
        pos(4) = 300;
        set(handles.figure1,'Position',pos);
    end
end

% New position of the axes
A = handles.axes;
cax = gca;
inx = find(A==cax);

nrows = handles.layout{1};  % get layout type
ncols = handles.layout{2};
layout(h,handles,nrows,ncols,A);

% Remember the old 'Units' properties and set 'Units' to pixels
old_units_axes = get(A(1),'Units');
set(A,'Units','pixels')

% New position of the pushbuttons
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
    
% Move legends
for x = 1:length(A)
    [L,legendpos] = legendcheck(handles.figure1,A(x));
    if ~isempty(L)
        old_axes_units = get(A(x),'Units');
        old_legend_units = get(L,'Units');
        set(A(x),'Units','points')
        set(L,'Units','points')
        pos = get(L,'Position');
        width = pos(3);
        hight = pos(4);
        lpos = legendposition(A(x),legendpos,width,hight);
        set(L,'Position',lpos,'Units',old_legend_units)
        set(A(x),'Units',old_axes_units)
    end
end

% Refresh 'handles' structure
handles = guidata(h);
A = handles.axes;
if ~isequal(gca,A(inx))
    error('1396')
end

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);
set(A,'Units',old_units_axes);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% --------------------------------------------------------------------
% LAYOUT CONSTRUCTION
% --------------------------------------------------------------------
function varargout = layout(h,handles,nrows,ncols,varargin)

% Input argument check
if nrows == 0 | ncols == 0
    handles.axes = [];
    guidata(h,handles)
    return
end
lenr = length(nrows);
lenc = length(ncols);
if lenr > 1 & lenc > 1
    error('Invalid input for layout construction: at least one of the inputs has to be scalar.')
end
if (lenr == 1 &~isequal(lenc,1) & ~isequal(nrows,lenc)) | (lenc == 1 &~isequal(lenr,1) & ~isequal(ncols,lenr))
    error('Invalid input for layout construction: dimension mismatch.')
end

% Store default properties
default_units = get(handles.figure1,'Units');
set(handles.figure1,'Units','pixels')
default_handlevisibility = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','on');

% Create axes
emp.name = '';
emp.settings = {};
emp.button = [];
emp.stringinput = '';
emp.keypress = '';
emp.opening_button = [];
if lenc > 1
    axno = sum(ncols);
elseif lenr > 1
    axno = sum(nrows);
else
    axno = nrows * ncols;
end
if isempty(varargin)
    for a = 1:axno
        A(a) = axes('Units','pixels');
    end
else
    A = varargin{1};
    set(A,'Units','pixels')
    lena = length(A);
    if lena > axno
        for a = axno+1:lena
            delete_legend2(A(a))
            delete_buttons2(A(a))
            delete(A(a))
        end
        A = A(1:axno);
        handles.guidisplay = handles.guidisplay(1:axno);
    elseif lena < axno
        for a = lena+1:axno
            A(a) = axes('Units','pixels');
            handles.guidisplay(a) = emp;
        end
    end
end

% Get GUI position
gui_pos = get(handles.figure1,'Position');
gui_left = gui_pos(1);
gui_bottom = gui_pos(2);
gui_width = gui_pos(3);
gui_hight = gui_pos(4);

% Set offset
offset = 40;

% Calculate axes positions and create axes
next = 0;
if lenc > 1
    horz = gui_hight - ((nrows + 1) * offset);
    axes_hight = horz / nrows;
    for r = 1:nrows
        vert(r) = gui_width - ((ncols(r) + 1) * offset);
        axes_width(r) = vert(r) / ncols(r);
        for c = 1:ncols(r)
            next = next + 1;
            axes_left = ((c - 1) * axes_width(r)) + (c * offset);
            axes_bottom = ((nrows - r) * axes_hight) + ((nrows - r + 1) * offset);
            axes_pos = [axes_left axes_bottom axes_width(r) axes_hight];
            set(A(next),'Position',axes_pos);
        end
    end
elseif lenr > 1
    vert = gui_width - ((ncols + 1) * offset);
    axes_width = vert / ncols;
    for c = 1:ncols
        horz(c) = gui_hight - ((nrows(c) + 1) * offset);
        axes_hight(c) = horz(c) / nrows(c);
        for r = 1:nrows(c)
            next = next + 1;
            axes_left = ((c - 1) * axes_width) + (c * offset);
            axes_bottom = ((nrows(c) - r) * axes_hight(c)) + ((nrows(c) - r + 1) * offset);
            axes_pos = [axes_left axes_bottom axes_width axes_hight(c)];
            set(A(next),'Position',axes_pos);
        end
    end
elseif lenr == 1 & lenc == 1
    vert = gui_width - ((ncols + 1) * offset);
    horz = gui_hight - ((nrows + 1) * offset);
    axes_width = vert / ncols;
    axes_hight = horz / nrows;
    for r = 1:nrows
        for c = 1:ncols
            next = next + 1;
            axes_left = ((c - 1) * axes_width) + (c * offset);
            axes_bottom = ((nrows - r) * axes_hight) + ((nrows - r + 1) * offset);
            axes_pos = [axes_left axes_bottom axes_width axes_hight];
            set(A(next),'Position',axes_pos);
        end
    end
end

% Set axes properties
set(A,'Box','on')
set(A,'FontSize',8)

% Buttondown function
bdf = 'b_ica_gui4(''Axes_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(A,'ButtonDownFcn',bdf)

% Store axes handles and layout type in 'handles' structure
handles.axes = A;
handles.layout = [{nrows} {ncols}];
guidata(h,handles)

% Uicontextmenus
cmenus(h,handles)

% Reset properties 
set(handles.figure1,'Units',default_units)
set(handles.figure1,'HandleVisibility',default_handlevisibility);
set(A,'Units','characters')



% --------------------------------------------------------------------
% UICONTEXTMENUS
% --------------------------------------------------------------------
function varargout = cmenus(h,handles)

% Store default properties
default_handlevisibility = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','on');

% Calculate special positions
[rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles);
handles = guidata(h);

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Get axes handles
A = handles.axes;

% Define the contextmenus
cmenu = cell(1,axno);
for i = 1:axno
    cmenu{i} = uicontextmenu;
end

% Define the context menu items
create_uimenus_common1(cmenu)     % common menu items
handles = guidata(h);
if lenr == 1 & lenc == 1
    create_uimenus_first(cmenu{1},'rowtype',nrows,ncols,handles.displen)      % menu item for first position
    create_uimenus_last(cmenu{end},'rowtype',nrows,ncols,handles.displen)      % menu item for last position
    create_uimenus_first(cmenu{1},'columntype',nrows,ncols,handles.displen)      % menu item for first position
    create_uimenus_last(cmenu{end},'columntype',nrows,ncols,handles.displen)      % menu item for last position
else
    create_uimenus_first(cmenu{1},type,nrows,ncols,handles.displen)      % menu item for first position
    create_uimenus_last(cmenu{end},type,nrows,ncols,handles.displen)      % menu item for last position    
end
create_uimenus_rowend({cmenu{rowend}},nrows,ncols,rowend,handles.displen)     % menu items for rowend position
create_uimenus_rowstart({cmenu{rowstart}},nrows,ncols,rowstart,handles.displen)     % menu items for rowstart position
create_uimenus_colend({cmenu{colend}},nrows,ncols,colend,handles.displen)     % menu items for colend position
create_uimenus_colstart({cmenu{colstart}},nrows,ncols,colstart,handles.displen)     % menu items for colstart position
create_uimenus_common2(cmenu)     % common menu items

% Associate the axes with the contexmenus
for i = 1:axno
    set(A(i),'UIContextMenu',cmenu{i})
end

% Reset properties 
set(handles.figure1,'HandleVisibility',default_handlevisibility);

% --------------------------------------------------------------------
function create_uimenus_common1(cmenu)

% Create context menu items
display_item_list(cmenu)

% --------------------------------------------------------------------
function create_uimenus_rowend(cmenu,nrows,ncols,rowend,displen)

% Define callbacks for context menu items
cb10 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''next_row'')';
cb11 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_lower_row'')';

% Define the context menu items
lenc = length(ncols);
if lenc == 1
    ncols(1:nrows) = ncols;
end
for i = 1:length(cmenu)
    for j = 1:lenc
        cumm(j) = sum(ncols(1:j));
    end
    fx = 0;
    fx = [fx find(cumm<rowend(i))];
    k = fx(end) + 1;
    if length(get(cmenu{i},'Children')) == displen
        uimenu(cmenu{i},'Label','Move to next row','Callback',cb10,'Separator','on');
    else
        uimenu(cmenu{i},'Label','Move to next row','Callback',cb10);
    end
    if ~isequal(ncols(k),1)
        uimenu(cmenu{i},'Label','Move to new lower row','Callback',cb11);
    end
end

% --------------------------------------------------------------------
function create_uimenus_rowstart(cmenu,nrows,ncols,rowstart,displen)

% Define callbacks for context menu items
cb12 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''previous_row'')';
cb13 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_upper_row'')';

% Define the context menu items
lenc = length(ncols);
if lenc == 1
    ncols(1:nrows) = ncols;
end
for i = 1:length(cmenu)
    for j = 1:lenc
        cumm(j) = sum(ncols(1:j));
    end
    fx = 0;
    fx = [fx find(cumm<rowstart(i))];
    k = fx(end) + 1;
    if length(get(cmenu{i},'Children')) == displen
        uimenu(cmenu{i},'Label','Move to previous row','Callback',cb12,'Separator','on');
    else
        uimenu(cmenu{i},'Label','Move to previous row','Callback',cb12);
    end
    if ~isequal(ncols(k),1)
        uimenu(cmenu{i},'Label','Move to new upper row','Callback',cb13);
    end
end

% --------------------------------------------------------------------
function create_uimenus_colend(cmenu,nrows,ncols,colend,displen)

% Define callbacks for context menu items
cb14 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''next_col'')';
cb15 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_right_col'')';

% Define the context menu items
lenr = length(nrows);
if lenr == 1
    nrows(1:ncols) = nrows;
end
for i = 1:length(cmenu)
    for j = 1:lenr
        cumm(j) = sum(nrows(1:j));
    end
    fx = 0;
    fx = [fx find(cumm<colend(i))];
    k = fx(end) + 1;
    if length(get(cmenu{i},'Children')) == displen
        uimenu(cmenu{i},'Label','Move to next column','Callback',cb14,'Separator','on');
    else
        uimenu(cmenu{i},'Label','Move to next column','Callback',cb14);
    end
    if ~isequal(nrows(k),1)
        uimenu(cmenu{i},'Label','Move to new right column','Callback',cb15);
    end
end

% --------------------------------------------------------------------
function create_uimenus_colstart(cmenu,nrows,ncols,colstart,displen)

% Define callbacks for context menu items
cb16 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''previous_col'')';
cb17 = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_left_col'')';

% Define the context menu items
lenr = length(nrows);
if lenr == 1
    nrows(1:ncols) = nrows;
end
for i = 1:length(cmenu)
    for j = 1:lenr
        cumm(j) = sum(nrows(1:j));
    end
    fx = 0;
    fx = [fx find(cumm<colstart(i))];
    k = fx(end) + 1;
    if length(get(cmenu{i},'Children')) == displen
        uimenu(cmenu{i},'Label','Move to previous column','Callback',cb16,'Separator','on');
    else
        uimenu(cmenu{i},'Label','Move to previous column','Callback',cb16);
    end
    if ~isequal(nrows(k),1)
        uimenu(cmenu{i},'Label','Move to new left column','Callback',cb17);
    end
end

% --------------------------------------------------------------------
function create_uimenus_common2(cmenu)

% Define callbacks for context menu items
cb18 = 'b_ica_gui4(''undock'',gcbo,guidata(gcbo))';

% Define the context menu items
for i = 1:length(cmenu)
    uimenu(cmenu{i},'Label','Undock','Callback',cb18,'Separator','on');
end

% --------------------------------------------------------------------
function create_uimenus_first(cmenu,type,nrows,ncols,displen)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_upper_row'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_left_col'')';

% Define the context menu items
if strcmp(type,'rowtype')
    if ~isequal(ncols(1),1)
        if length(get(cmenu,'Children')) == displen
            uimenu(cmenu,'Label','Move to new upper row','Callback',cb_a,'Separator','on');
        else
            uimenu(cmenu,'Label','Move to new upper row','Callback',cb_a);
        end
    end
elseif strcmp(type,'columntype')
    if ~isequal(nrows(1),1)
        if length(get(cmenu,'Children')) == displen
            uimenu(cmenu,'Label','Move to new left column','Callback',cb_b,'Separator','on');
        else
            uimenu(cmenu,'Label','Move to new left column','Callback',cb_b);
        end
    end
end

% --------------------------------------------------------------------
function create_uimenus_last(cmenu,type,nrows,ncols,displen)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_lower_row'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_right_col'')';

% Define the context menu items
if strcmp(type,'rowtype')
    if ~isequal(ncols(end),1)
        if length(get(cmenu,'Children')) == displen
            uimenu(cmenu,'Label','Move to new lower row','Callback',cb_a,'Separator','on');
        else
            uimenu(cmenu,'Label','Move to new lower row','Callback',cb_a);
        end
    end
elseif strcmp(type,'columntype')
    if ~isequal(nrows(end),1)
        if length(get(cmenu,'Children')) == displen
            uimenu(cmenu,'Label','Move to new right column','Callback',cb_b,'Separator','on');
        else
            uimenu(cmenu,'Label','Move to new right column','Callback',cb_b)
        end
    end
end

% --------------------------------------------------------------------
function varargout = figure_cmenus(h,handles,ax)

% Calculate special positions
[rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles);
handles = guidata(h);
    
% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Define the contextmenu
cmenu = uicontextmenu('Parent',get(ax,'Parent'));

% Associate the axes with the contexmenu
set(ax,'UIContextMenu',cmenu)
ach = get(ax,'Children');
set(ach,'UIContextMenu',cmenu)

% Delete buttondown functions
set(ach,'ButtonDownFcn','');
set(ax,'ButtonDownFcn','');

% Define callbacks for context menu items
posno = axno + 1;
cb = cell(1,posno+1);
cb{1} = 'closereq';
for i = 1:posno
    cb{i+1} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(i) ','''')'];
end

% Define the context menu items
flag{1} = 'st';
flag{2} = 'nd';
flag{3} = 'rd';
if posno > 3
    for i = 4:posno
        flag{i} = 'th';
    end
end
item = cell(1,posno+1);
item{1} = uimenu(cmenu,'Label','Close','Callback',cb{1});
for i = 1:posno
    item{i+1} = uimenu(cmenu,'Label',['Dock to ' num2str(i) flag{i} ' position'],...
        'Callback',cb{i+1});
end

% Create submenus for accurate docking
create_figcmenus_first(item{2},type)      % submenu item for first position
create_figcmenus_last(item{end},type,nrows,ncols)      % submenu item for last position
if strcmp(type,'rowtype')
    create_figcmenus_rowstart({item{rowstart+1}},nrows,ncols,rowstart)     % submenu items for rowstart position
elseif strcmp(type,'columntype')
    create_figcmenus_colstart({item{colstart+1}},nrows,ncols,colstart)     % submenu items for colstart position
end

% --------------------------------------------------------------------
function create_figcmenus_rowstart(item,nrows,ncols,rowstart)

% Define callbacks for context menu items
leni = length(item);
cb1 = cell(1,leni);
cb2 = cell(1,leni);
cb3 = cell(1,leni);
for i = 1:length(item)
    cb1{i} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(rowstart(i)) ',''upper_row'')'];
    cb2{i} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(rowstart(i)) ',''new_row'')'];
    cb3{i} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(rowstart(i)) ',''lower_row'')'];
end

% Define the context menu items
lenc = length(ncols);
if lenc == 1
    ncols(1:nrows) = ncols;
end
for i = 1:length(item)
    uimenu(item{i},'Label','Upper row','Callback',cb1{i});
    uimenu(item{i},'Label','New row','Callback',cb2{i});
    uimenu(item{i},'Label','Lower row','Callback',cb3{i});
    set(item{i},'Callback','')      % remove original contextmenu callback
end

% --------------------------------------------------------------------
function create_figcmenus_colstart(item,nrows,ncols,colstart)

% Define callbacks for context menu items
leni = length(item);
cb1 = cell(1,leni);
cb2 = cell(1,leni);
cb3 = cell(1,leni);
for i = 1:length(item)
    cb1{i} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(colstart(i)) ',''left_col'')'];
    cb2{i} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(colstart(i)) ',''new_col'')'];
    cb3{i} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(colstart(i)) ',''right_col'')'];
end

% Define the context menu items
lenr = length(nrows);
if lenr == 1
    nrows(1:ncols) = nrows;
end
for i = 1:length(item)
    uimenu(item{i},'Label','Left column','Callback',cb1{i});
    uimenu(item{i},'Label','New column','Callback',cb2{i});
    uimenu(item{i},'Label','Right column','Callback',cb3{i});
    set(item{i},'Callback','')      % remove original contextmenu callback
end

% --------------------------------------------------------------------
function create_figcmenus_first(item,type)

% Define callbacks for context menu items
cb_a1 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),1,''new_first_row'')'];
cb_a2 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),1,''first_row'')'];
cb_b1 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),1,''new_first_col'')'];
cb_b2 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),1,''first_col'')'];

% Define the context menu items
if strcmp(type,'rowtype')
    uimenu(item,'Label','New row','Callback',cb_a1);
    uimenu(item,'Label','First row','Callback',cb_a2);
elseif strcmp(type,'columntype')
    uimenu(item,'Label','New column','Callback',cb_b1);
    uimenu(item,'Label','First column','Callback',cb_b2);
end

% Remove original contextmenu callback
set(item,'Callback','')

% --------------------------------------------------------------------
function create_figcmenus_last(item,type,nrows,ncols)

% Define callbacks for context menu items
if length(nrows) == 1 & length(ncols) == 1
    lastpos = nrows * ncols + 1;
else
    if strcmp(type,'rowtype')
        lastpos = sum(ncols) + 1;
    elseif strcmp(type,'columntype')
        lastpos = sum(nrows) + 1;
    end
end
cb_a1 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(lastpos) ',''new_last_row'')'];
cb_a2 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(lastpos) ',''last_row'')'];
cb_b1 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(lastpos) ',''new_last_col'')'];
cb_b2 = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(lastpos) ',''last_col'')'];

% Define the context menu items
if strcmp(type,'rowtype')
    uimenu(item,'Label','New row','Callback',cb_a1);
    uimenu(item,'Label','Last row','Callback',cb_a2);
elseif strcmp(type,'columntype')
    uimenu(item,'Label','New column','Callback',cb_b1);
    uimenu(item,'Label','Last column','Callback',cb_b2);
end

% Remove original contextmenu callback
set(item,'Callback','')



% --------------------------------------------------------------------
% UIMENUS
% --------------------------------------------------------------------
function varargout = menus(h,handles)

% Store default properties
default_handlevisibility = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','on');

% Calculate special positions
[rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles);
handles = guidata(h);

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Get axes handles
cax = gca;
A = handles.axes;
inx = find(A==cax);

% Define the menus
vmenu = handles.view_menu;
dmenu = handles.display_menu;

% Create menu items for 'Display' menu
display_item_list(dmenu)

% Define callbacks for 'View' menu items
cb10 = 'b_ica_layouteditor';
cb11 = 'b_ica_gui4(''typeswitch'',gcbo,guidata(gcbo))';
cb12 = 'b_ica_gui4(''undock'',gcbo,guidata(gcbo))';

% Define the 'View' menu items
item10 = uimenu(vmenu,'Label','Layout','Callback',cb10);
if (isequal(lenc,1) & isequal(lenr,1)) | (lenr > 1 & b_isconstant(nrows)) | (lenc > 1 & b_isconstant(ncols))
    str = 'on';
else
    str = 'off';
end
item11 = uimenu(vmenu,'Label','Switch type','Callback',cb11,'Enable',str,'Accelerator','y');
item12 = uimenu(vmenu,'Label','Undock', 'Callback',cb12,'Separator','on','Accelerator','u');
if isequal(inx,1)
    if lenr == 1 & lenc == 1
        create_submenu_first(vmenu,'rowtype',nrows,ncols)      % menu item for first position
        create_submenu_first(vmenu,'columntype',nrows,ncols)      % menu item for first position
    else
        create_submenu_first(vmenu,type,nrows,ncols);      % menu item for first position
    end
end
if isequal(inx,length(A))
    if lenr == 1 & lenc == 1
        create_submenu_last(vmenu,'rowtype',nrows,ncols)      % menu item for last position
        create_submenu_last(vmenu,'columntype',nrows,ncols)      % menu item for last position
    else
        create_submenu_last(vmenu,type,nrows,ncols);      % menu item for last position
    end
end
if ~isempty(find(([0.5 rowend]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
    create_submenu_rowend(vmenu,nrows,ncols,inx)     % menu items for rowend position
end
if ~isempty(find(([0.5 rowstart]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
    create_submenu_rowstart(vmenu,nrows,ncols,inx)     % menu items for rowstart position
end
if ~isempty(find(([0.5 colend]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
    create_submenu_colend(vmenu,nrows,ncols,inx)     % menu items for colend position
end
if ~isempty(find(([0.5 colstart]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
    create_submenu_colstart(vmenu,nrows,ncols,inx)     % menu items for colstart position
end

% Create shortcuts for submenus in 'File' menu
set(handles.open_submenu,'Accelerator','o')
set(handles.next_submenu,'Accelerator','n')
set(handles.previous_submenu,'Accelerator','b')

% Reset properties 
set(handles.figure1,'HandleVisibility',default_handlevisibility);

% --------------------------------------------------------------------
function create_submenu_rowend(vmenu,nrows,ncols,inx)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''next_row'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_lower_row'')';

% Define the context menu items
lenc = length(ncols);
if lenc == 1
    ncols(1:nrows) = ncols;
end
for j = 1:lenc
    cumm(j) = sum(ncols(1:j));
end
fx = 0;
fx = [fx find(cumm<inx)];
k = fx(end) + 1;
if length(get(vmenu,'Children')) == 1
    uimenu(vmenu,'Label','Move to next row','Callback',cb_a,'Separator','on','Accelerator','r');
else
    uimenu(vmenu,'Label','Move to next row','Callback',cb_a,'Accelerator','r');
end
if ~isequal(ncols(k),1)
    uimenu(vmenu,'Label','Move to new lower row','Callback',cb_b,'Accelerator','d');
end

% --------------------------------------------------------------------
function create_submenu_rowstart(vmenu,nrows,ncols,inx)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''previous_row'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_upper_row'')';

% Define the context menu items
lenc = length(ncols);
if lenc == 1
    ncols(1:nrows) = ncols;
end
for j = 1:lenc
    cumm(j) = sum(ncols(1:j));
end
fx = 0;
fx = [fx find(cumm<inx)];
k = fx(end) + 1;
if length(get(vmenu,'Children')) == 1
    uimenu(vmenu,'Label','Move to previous row','Callback',cb_a,'Separator','on','Accelerator','t');
else
    uimenu(vmenu,'Label','Move to previous row','Callback',cb_a,'Accelerator','t');
end
if ~isequal(ncols(k),1)
    uimenu(vmenu,'Label','Move to new upper row','Callback',cb_b,'Accelerator','f');
end

% --------------------------------------------------------------------
function create_submenu_colend(vmenu,nrows,ncols,inx)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''next_col'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_right_col'')';

% Define the context menu items
lenr = length(nrows);
if lenr == 1
    nrows(1:ncols) = nrows;
end
for j = 1:lenr
    cumm(j) = sum(nrows(1:j));
end
fx = 0;
fx = [fx find(cumm<inx)];
k = fx(end) + 1;
if length(get(vmenu,'Children')) == 1
    uimenu(vmenu,'Label','Move to next column','Callback',cb_a,'Separator','on','Accelerator','w');
else
    uimenu(vmenu,'Label','Move to next column','Callback',cb_a,'Accelerator','w');
end
if ~isequal(nrows(k),1)
    uimenu(vmenu,'Label','Move to new right column','Callback',cb_b,'Accelerator','a');
end

% --------------------------------------------------------------------
function create_submenu_colstart(vmenu,nrows,ncols,inx)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''previous_col'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_left_col'')';

% Define the context menu items
lenr = length(nrows);
if lenr == 1
    nrows(1:ncols) = nrows;
end
for j = 1:lenr
    cumm(j) = sum(nrows(1:j));
end
fx = 0;
fx = [fx find(cumm<inx)];
k = fx(end) + 1;
if length(get(vmenu,'Children')) == 1
    uimenu(vmenu,'Label','Move to previous column','Callback',cb_a,'Separator','on','Accelerator','e');
else
    uimenu(vmenu,'Label','Move to previous column','Callback',cb_a,'Accelerator','e');
end
if ~isequal(nrows(k),1)
    uimenu(vmenu,'Label','Move to new left column','Callback',cb_b,'Accelerator','s');
end

% --------------------------------------------------------------------
function create_submenu_first(vmenu,type,nrows,ncols)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_upper_row'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_left_col'')';

% Define the context menu items
if strcmp(type,'rowtype')
    if ~isequal(ncols(1),1)
        if length(get(vmenu,'Children')) == 1
            uimenu(vmenu,'Label','Move to new upper row','Callback',cb_a,'Separator','on','Accelerator','f');
        else
            uimenu(vmenu,'Label','Move to new upper row','Callback',cb_a,'Accelerator','f');
        end
    end
elseif strcmp(type,'columntype')
    if ~isequal(nrows(1),1)
        if length(get(vmenu,'Children')) == 1
            uimenu(vmenu,'Label','Move to new left column','Callback',cb_b,'Separator','on','Accelerator','s');
        else
            uimenu(vmenu,'Label','Move to new left column','Callback',cb_b,'Accelerator','s');
        end
    end
end

% --------------------------------------------------------------------
function create_submenu_last(vmenu,type,nrows,ncols)

% Define callbacks for context menu items
cb_a = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_lower_row'')';
cb_b = 'b_ica_gui4(''move_axes'',gcbo,guidata(gcbo),''new_right_col'')';

% Define the context menu items
if strcmp(type,'rowtype')
    if ~isequal(ncols(end),1)
        if length(get(vmenu,'Children')) == 1
            uimenu(vmenu,'Label','Move to new lower row','Callback',cb_a,'Separator','on','Accelerator','d');
        else
            uimenu(vmenu,'Label','Move to new lower row','Callback',cb_a,'Accelerator','d');
        end
    end
elseif strcmp(type,'columntype')
    if ~isequal(nrows(end),1)
        if length(get(vmenu,'Children')) == 1
            uimenu(vmenu,'Label','Move to new right column','Callback',cb_b,'Separator','on','Accelerator','a');
        else
            uimenu(vmenu,'Label','Move to new right column','Callback',cb_b,'Accelerator','a')
        end
    end
end

% --------------------------------------------------------------------
function varargout = figure_menus(h,handles)

% Calculate the special positions
[rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles);
handles = guidata(h);

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
lenr = length(nrows);
lenc = length(ncols);

% Define the menus
vmenu = uimenu('Label','View');
dmenu = uimenu('Label','Display');
tmenu = uimenu('Label','Tools');

% Create 'Display' menu items
display_item_list(dmenu)

% Define callbacks for 'View' menu items
posno = axno + 1;
cb = cell(1,posno+1);
cb{1} = 'closereq';
for i = 1:posno
    cb{i+1} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(i) ','''')'];
end

% Define the menu items for 'View' menu
flag{1} = 'st';
flag{2} = 'nd';
flag{3} = 'rd';
if posno > 3
    for i = 4:posno
        flag{i} = 'th';
    end
end

item = cell(1,posno+1);
item{1} = uimenu(vmenu,'Label','Close','Callback',cb{1},'Accelerator','q');
for i = 1:min(posno,9)    % only the 1-9 position has a shortcut
    item{i+1} = uimenu(vmenu,'Label',['Dock to ' num2str(i) flag{i} ' position'],...
        'Callback',cb{i+1},'Accelerator',[num2str(i)]);
end
if posno > 9 
    for i = 10:posno
        item{i+1} = uimenu(vmenu,'Label',['Dock to ' num2str(i) flag{i} ' position'],...
            'Callback',cb{i+1});
    end
end

% Define callbacks for 'Tools' menu items
zs = 'b_ica_gui4(''figure_zoom'');';

% Define the menu items for 'Tools' menu
item1 = uimenu(tmenu,'Label','Zoom','Callback',zs);

% --------------------------------------------------------------------
function display_item_list(dmenu)

% Convert input to cell
if ~iscell(dmenu)
    temp{1} = dmenu;
    dmenu = temp;
end

% Define callbacks for menu items
dl = 'b_ica_gui4(''delete_legend'',gcbo,guidata(gcbo));';
db = 'b_ica_gui4(''delete_buttons'',gcbo,guidata(gcbo));';
dd = [dl db];
cb1 = [dd 'b_ica_gui4(''variance_plot'',gcbo,guidata(gcbo))'];
cb2 = [dd 'b_ica_gui4(''structural_elements'',gcbo,guidata(gcbo),1)'];
cb3 = [dd 'b_ica_gui4(''return_plots'',gcbo,guidata(gcbo),1)'];
cb4 = [dd 'b_ica_gui4(''burst_parameters'',gcbo,guidata(gcbo),1)'];
cb4p1 = [dd 'b_ica_gui4(''instantenous_frequency'',gcbo,guidata(gcbo),1)'];
% cb5 = [dd 'b_ica_gui4(''lomb_periodogram'',gcbo,guidata(gcbo),1)'];
cb6a = [dd 'b_ica_gui4(''wavelet'',gcbo,guidata(gcbo),''power'')'];
cb6b = [dd 'b_ica_gui4(''wavelet'',gcbo,guidata(gcbo),''phase'')'];
cb7 = [dd 'b_ica_gui4(''wavelet_average'',gcbo,guidata(gcbo))'];
cb8 = [dd 'b_ica_gui4(''wavelet_average_fft'',gcbo,guidata(gcbo),1)'];
cb9a = [dd 'b_ica_gui4(''burstwave_matrices'',gcbo,guidata(gcbo),1,''power'')'];
cb9b = [dd 'b_ica_gui4(''burstwave_matrices'',gcbo,guidata(gcbo),1,''phase'')'];
cb10a = [dd 'b_ica_gui4(''scale_cv_functions'',gcbo,guidata(gcbo),''power'')'];
cb10b = [dd 'b_ica_gui4(''scale_cv_functions'',gcbo,guidata(gcbo),''phase'')'];
cb11a = [dd 'b_ica_gui4(''spike_triggered_wavelet'',gcbo,guidata(gcbo),1,''power'')'];
cb11b = [dd 'b_ica_gui4(''spike_triggered_wavelet'',gcbo,guidata(gcbo),1,''phase'')'];
cb11c = [dd 'b_ica_gui4(''spike_triggered_wavelet_average'',gcbo,guidata(gcbo),''power'')'];
cb11d = [dd 'b_ica_gui4(''spike_triggered_wavelet_average'',gcbo,guidata(gcbo),''phase'')'];
cb12a = [dd 'b_ica_gui4(''spstw'',gcbo,guidata(gcbo),1,''power'')'];
cb12b = [dd 'b_ica_gui4(''spstw'',gcbo,guidata(gcbo),1,''phase'')'];
cb12c = [dd 'b_ica_gui4(''spstw_average'',gcbo,guidata(gcbo),1,''power'')'];
cb12d = [dd 'b_ica_gui4(''spstw_average'',gcbo,guidata(gcbo),1,''phase'')'];
cb13 = [dd 'b_ica_gui4(''thresholding'',gcbo,guidata(gcbo),1)'];

% Define the menu items
for i = 1:length(dmenu)
    item1 = uimenu(dmenu{i},'Label','Variance plot','Callback',cb1);
    item2 = uimenu(dmenu{i},'Label','Structural elements','Callback',cb2);
    item3 = uimenu(dmenu{i},'Label','Return plots','Callback',cb3);
    item4 = uimenu(dmenu{i},'Label','Burst parameters','Callback',cb4);
    item4p1 = uimenu(dmenu{i},'Label','Instantenous frequency','Callback',cb4p1);
%     item5 = uimenu(dmenu{i},'Label','Lomb periodogram','Callback',cb5);
    item6 = uimenu(dmenu{i},'Label','Wavelet');
    item6a = uimenu(item6,'Label','Power','Callback',cb6a);
    item6b = uimenu(item6,'Label','Phase','Callback',cb6b);
    item7 = uimenu(dmenu{i},'Label','Wavelet average','Callback',cb7);
    item8 = uimenu(dmenu{i},'Label','Wavelet average FFT','Callback',cb8);
    item9 = uimenu(dmenu{i},'Label','Burstwave matrices');
    item9a = uimenu(item9,'Label','Power','Callback',cb9a);
    item9b = uimenu(item9,'Label','Phase','Callback',cb9b);
    item10 = uimenu(dmenu{i},'Label','Scale-CV functions');
    item10a = uimenu(item10,'Label','Power','Callback',cb10a);
    item10b = uimenu(item10,'Label','Phase','Callback',cb10b);
    item11 = uimenu(dmenu{i},'Label','Spike triggered wavelet');
    item11a = uimenu(item11,'Label','Power','Callback',cb11a);
    item11b = uimenu(item11,'Label','Phase','Callback',cb11b);
    item11c = uimenu(item11,'Label','Power, averaged','Callback',cb11c);
    item11d = uimenu(item11,'Label','Phase, averaged','Callback',cb11d);
    item12 = uimenu(dmenu{i},'Label','Single point spike triggered wavelet');
    item12a = uimenu(item12,'Label','Power','Callback',cb12a);
    item12b = uimenu(item12,'Label','Phase','Callback',cb12b);
    item12c = uimenu(item12,'Label','Power, averaged','Callback',cb12c);
    item12d = uimenu(item12,'Label','Phase, averaged','Callback',cb12d);
    item13 = uimenu(dmenu{i},'Label','Thresholding','Callback',cb13);
end

% Length of display list
global ICA_GUI_HANDLE
h = ICA_GUI_HANDLE;
handles = guidata(h);
handles.displen = 13;
guidata(h,handles)



% --------------------------------------------------------------------
% DOCK & UNDOCK
% --------------------------------------------------------------------
function undock(h,handles)

% Get axes handles
cax = gca;
A = handles.axes;
inx = find(A==cax);
ux = A(inx);

% Change 'h' to gui handle
h = handles.figure1;

% New layout
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);
if strcmp(type,'rowtype')
    if lenc == 1
        ncols(1:nrows) = ncols;    % transform types different from one row
        lenc = nrows;
    end
    if lenc > 1
        for i = 1:lenc
            cumm(i) = sum(ncols(1:i));
        end
        fx = 0;
        fx = [fx find(cumm<inx)];
        nd = fx(end) + 1;
        ncols(nd) = ncols(nd) - 1;
        if ncols(nd) == 0
            nrows = nrows - 1;
            ncols(nd) = [];
        end
    elseif lenc == 1
        ncols = ncols - 1;    % one row case
    end
elseif strcmp(type,'columntype')
    if lenr == 1
        nrows(1:ncols) = nrows;     % transform types different from one column
        lenr = ncols;
    end
    if lenr > 1
        for i = 1:lenr
            cumm(i) = sum(nrows(1:i));
        end
        fx = 0;
        fx = [fx find(cumm<inx)];
        nd = fx(end) + 1;
        nrows(nd) = nrows(nd) - 1;
        if nrows(nd) == 0
            ncols = ncols - 1;
            nrows(nd) = [];
        end
    elseif lenr == 1
        nrows = nrows - 1;     % one column case
    end
end
A(inx) = [];
layout(h,handles,nrows,ncols,A)

% Refresh 'handles' structure
handles = guidata(h);

% New figure
H = figure('MenuBar','none','ResizeFcn','b_ica_gui4(''resize_undocked'')','CloseRequestFcn',...
    'b_ica_gui4(''closereq_undocked'')');
set(ux,'Units',get(H,'DefaultAxesUnits'));
set(ux,'Parent',H,'Position',get(H,'DefaultAxesPosition'),'FontSize',10,'XColor','black','YColor','Black',...
    'ButtonDownFcn',[])
handles.figure(end+1) = H;
guidata(h,handles)

% Move legend
[L,legendpos] = legendcheck(handles.figure1,ux);
if ~isempty(L)
    old_axes_units = get(ux,'Units');
    old_legend_units = get(L,'Units');
    set(ux,'Units','points')
    set(L,'Units','points')
    pos = get(L,'Position');
    width = pos(3);
    hight = pos(4);
    lpos = legendposition(ux,legendpos,width,hight);
    set(L,'Parent',H,'Position',lpos,'Units',old_legend_units)
    set(ux,'Units',old_axes_units)
end

% Move buttons
button_handles = findobj(get(handles.figure1,'Children'),'Style','pushbutton');
button_handles = flipud(button_handles);
if ~isempty(button_handles)     % remember old 'Units'
    old_units_buttons = get(button_handles(1),'Units');
end
old_axes_units = get(A(1),'Units');
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
    default_units = get(ax,'Units');
    set(ax,'Units','pixels')
    lenb = usd.button_number;
    ax_pos = get(ax,'Position');
    width = ax_pos(3) / lenb;
    left = ax_pos(1) + (next - 1) * width;
    buttom = ax_pos(2) + ax_pos(4);
    hight = 20;
    buttons_pos = [left buttom width hight];
    if isequal(ax,cax)
        set(button_handles(i),'Parent',H,'Position',buttons_pos);
    else
        set(button_handles(i),'Position',buttons_pos);
    end
    set(ax,'Units',default_units)
    next = next + 1;
end

if ~isempty(button_handles)
    set(button_handles,'Units',old_units_buttons);   % reset 'Units'
end
set(A,'Units',old_axes_units)

% Figure menus & context menus
figure_menus(h,handles)
figure_cmenus(h,handles,cax)

% Refresh 'View' GUI menu
menu_refresher(h,handles);

% Refresh 'display' fields
emp.name = '';
emp.settings = {};
emp.button = [];
emp.stringinput = '';
emp.keypress = '';
emp.opening_button = [];
if length(handles.guidisplay) < inx
    handles.figdisplay(end+1) = emp;
else
    handles.figdisplay(end+1) = handles.guidisplay(inx);
    handles.guidisplay(inx) = [];
end
guidata(h,handles)

% Set keypress function
set_keypressfcn(h,handles)

% Refresh figure menus & context menus
for i = 1:length(handles.figure)
    refresh_figure_menus(h,handles,handles.figure(i))
    ax = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
    figure_cmenus(h,handles,ax)
end

% --------------------------------------------------------------------
function dock(h,handles,position,str)

% Get figure handle
H = gcf;

% Get 'handles' structure
global ICA_GUI_HANDLE
h = ICA_GUI_HANDLE;
handles = guidata(h);

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Get axes handles
cax = gca;
A = handles.axes;
for i = 1:length(handles.figure)
    B(i) = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
end
inx = find(B==cax);

% Calculate 'current row/column'
if strcmp(type,'rowtype')
    if lenc == 1
        ncols(1:nrows) = ncols;    % transform types different from one row
        lenc = nrows;
    end
    if lenc > 1
        for i = 1:lenc
            cumm(i) = sum(ncols(1:i));
        end
        fx = 0;
        fx = [fx find(cumm<position)];
        currow = fx(end) + 1;
    elseif lenc == 1
        currow = 1;     % one row case
    end
elseif strcmp(type,'columntype')
    if lenr == 1
        nrows(1:ncols) = nrows;     % transform types different from one column
        lenr = ncols;
    end
    if lenr > 1
        for i = 1:lenr
            cumm(i) = sum(nrows(1:i));
        end
        fx = 0;
        fx = [fx find(cumm<position)];
        curcol = fx(end) + 1;
    elseif lenr == 1
        curcol = 1;     % one column case
    end
end

% New layout
switch str
case ''
    if strcmp(type,'rowtype')
        if ~isequal(position,sum(ncols)+1)
            ncols(currow) = ncols(currow) + 1;
        else
            ncols(end) = ncols(end) + 1;
        end
    elseif strcmp(type,'columntype')
        if ~isequal(position,sum(nrows)+1)
            nrows(curcol) = nrows(curcol) + 1;
        else
            nrows(end) = nrows(end) + 1;
        end
    end
case 'lower_row'
    ncols(currow) = ncols(currow) + 1;
case 'new_row'
    ncols = [ncols(1:currow-1) 1 ncols(currow:end)];
    nrows = nrows + 1;
case 'upper_row'
    ncols(currow-1) = ncols(currow-1) + 1;
case 'first_row'
    ncols(1) = ncols(1) + 1;
case 'new_first_row'
    ncols = [1 ncols];
    nrows = nrows + 1;
case 'last_row'
    ncols(end) = ncols(end) + 1;
case 'new_last_row'
    ncols = [ncols 1];
    nrows = nrows + 1;
case 'right_col'
    nrows(curcol) = nrows(curcol) + 1;
case 'new_col'
    nrows = [nrows(1:curcol-1) 1 nrows(curcol:end)];
    ncols = ncols + 1;
case 'left_col'
    nrows(curcol-1) = nrows(curcol-1) + 1;
case 'first_col'
    nrows(1) = nrows(1) + 1;
case 'new_first_col'
    nrows = [1 nrows];
    ncols = ncols + 1;
case 'last_col'
    nrows(end) = nrows(end) + 1;
case 'new_last_col'
    nrows = [nrows 1];
    ncols = ncols + 1;
end
figure(h)
bdf = 'b_ica_gui4(''Axes_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(cax,'Parent',handles.figure1,'FontSize',8,'XColor','red','YColor','red','ButtonDownFcn',bdf)
A = [A(1:position-1) cax A(position:end)];
layout(h,handles,nrows,ncols,A)

% Refresh 'handles' structure
handles = guidata(h);

% Move legend
[L,legendpos] = legendcheck(H,cax);
if ~isempty(L)
    old_axes_units = get(cax,'Units');
    old_legend_units = get(L,'Units');
    set(cax,'Units','points')
    set(L,'Units','points')
    pos = get(L,'Position');
    width = pos(3);
    hight = pos(4);
    lpos = legendposition(cax,legendpos,width,hight);
    set(L,'Parent',h,'Position',lpos,'Units',old_legend_units)
    set(cax,'Units',old_axes_units)
end

% Move buttons1 - buttons on undocked figure
button_handles = findobj(get(H,'Children'),'Style','pushbutton');
button_handles = flipud(button_handles);
if ~isempty(button_handles)     % remember old 'Units'
    old_units_buttons = get(button_handles(1),'Units');
end
old_axes_units = get(A(1),'Units');
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
    default_units = get(ax,'Units');
    set(ax,'Units','pixels')
    lenb = usd.button_number;
    ax_pos = get(ax,'Position');
    width = ax_pos(3) / lenb;
    left = ax_pos(1) + (next - 1) * width;
    buttom = ax_pos(2) + ax_pos(4);
    hight = 20;
    buttons_pos = [left buttom width hight];
    set(button_handles(i),'Parent',handles.figure1,'Position',buttons_pos);
    set(ax,'Units',default_units)
    next = next + 1;
end

if ~isempty(button_handles)
    set(button_handles,'Units',old_units_buttons);   % reset 'Units'
end
set(A,'Units',old_axes_units)

% Move buttons2 - buttons on GUI
button_handles = findobj(get(handles.figure1,'Children'),'Style','pushbutton');
button_handles = flipud(button_handles);
if ~isempty(button_handles)     % remember old 'Units'
    old_units_buttons = get(button_handles(1),'Units');
end
old_axes_units = get(A(1),'Units');
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

% Refresh 'display' fields
handles.guidisplay = [handles.guidisplay(1:position-1) handles.figdisplay(inx) handles.guidisplay(position:end)];
handles.figdisplay(inx) = [];
handles.figure(inx) = [];
guidata(h,handles)

% Delete figure
delete(H)

% Refresh 'View' menu
menu_refresher(h,handles);

% Refresh figure menus & context menus
for i = 1:length(handles.figure)
    refresh_figure_menus(h,handles,handles.figure(i))
    ax = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
    figure_cmenus(h,handles,ax)
end



% --------------------------------------------------------------------
% LAYOUT REORGANISATION
% --------------------------------------------------------------------
function move_axes(h,handles,str)

% Gain 'h'
h = handles.figure1;

% Get layout type - part 1
type = handles.type;

% Automatic typeswitch if needed
lgc1 = strcmp(type,'rowtype') & ~isempty(find(strcmp(str,{'next_col','previous_col','new_left_col','new_right_col'})));
lgc2 = strcmp(type,'columntype') & ~isempty(find(strcmp(str,{'next_row','previous_row','new_upper_row','new_lower_row'})));
lgc3 = lgc1 | lgc2;
if lgc3
    typeswitch(h,handles)
    handles = guidata(h);     % refresh 'handles' structure
    type = handles.type;      % refresh layout type
end

% Get layout type - part 2
nrows = handles.layout{1};
ncols = handles.layout{2};
lenr = length(nrows);
lenc = length(ncols);

% Get axes handles
cax = gca;
A = handles.axes;
inx = find(A==cax);

% Calculate 'current row/column'
if strcmp(type,'rowtype')
    if lenc == 1
        ncols(1:nrows) = ncols;    % transform types different from one row
        lenc = nrows;
    end
    if lenc > 1
        for i = 1:lenc
            cumm(i) = sum(ncols(1:i));
        end
        fx = 0;
        fx = [fx find(cumm<inx)];
        currow = fx(end) + 1;
    elseif lenc == 1
        currow = 1;     % one row case
    end
elseif strcmp(type,'columntype')
    if lenr == 1
        nrows(1:ncols) = nrows;     % transform types different from one column
        lenr = ncols;
    end
    if lenr > 1
        for i = 1:lenr
            cumm(i) = sum(nrows(1:i));
        end
        fx = 0;
        fx = [fx find(cumm<inx)];
        curcol = fx(end) + 1;
    elseif lenr == 1
        curcol = 1;     % one column case
    end
end

% New layout
switch str
case {'next_row','previous_row','new_upper_row','new_lower_row'}
    ncols(currow) = ncols(currow) - 1;
case {'next_col','previous_col','new_left_col','new_right_col'}
    nrows(curcol) = nrows(curcol) - 1;
end

switch str
case 'next_row'
    ncols(currow+1) = ncols(currow+1) + 1;
case 'previous_row'
    ncols(currow-1) = ncols(currow-1) + 1;
case 'next_col'
    nrows(curcol+1) = nrows(curcol+1) + 1;
case 'previous_col'
    nrows(curcol-1) = nrows(curcol-1) + 1;
case 'new_upper_row'
    ncols = [ncols(1:currow-1) 1 ncols(currow:end)];
    nrows = nrows + 1;
case 'new_lower_row'
    ncols = [ncols(1:currow) 1 ncols(currow+1:end)];
    nrows = nrows + 1;
case 'new_left_col'
    nrows = [nrows(1:curcol-1) 1 nrows(curcol:end)];
    ncols = ncols + 1;
case 'new_right_col'
    nrows = [nrows(1:curcol) 1 nrows(curcol+1:end)];
    ncols = ncols + 1;
end

switch str
case {'next_row','previous_row','new_upper_row','new_lower_row'}
    if ncols(currow) == 0
        ncols(currow) = [];
        nrows = nrows - 1;
    end
case {'next_col','previous_col','new_left_col','new_right_col'}
    if nrows(curcol) == 0
        nrows(curcol) = [];
        ncols = ncols - 1;
    end
end

layout(h,handles,nrows,ncols,A)

% Move legends
for x = 1:length(A)
    [L,legendpos] = legendcheck(handles.figure1,A(x));
    if ~isempty(L)
        old_axes_units = get(A(x),'Units');
        old_legend_units = get(L,'Units');
        set(A(x),'Units','points')
        set(L,'Units','points')
        pos = get(L,'Position');
        width = pos(3);
        hight = pos(4);
        lpos = legendposition(A(x),legendpos,width,hight);
        set(L,'Position',lpos,'Units',old_legend_units)
        set(A(x),'Units',old_axes_units)
    end
end

% Move buttons
button_handles = findobj(get(handles.figure1,'Children'),'Style','pushbutton');
button_handles = flipud(button_handles);
if ~isempty(button_handles)     % remember old 'Units'
    old_units_buttons = get(button_handles(1),'Units');
end
old_axes_units = get(A(1),'Units');
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

% Refresh 'handles' structure
handles = guidata(h);
A = handles.axes;
if ~isequal(gca,A(inx))
    error('1396')
end

% Refresh 'View' menu
menu_refresher(h,handles);

% Refresh figure menus & context menus
for i = 1:length(handles.figure)
    refresh_figure_menus(h,handles,handles.figure(i))
    ax = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
    figure_cmenus(h,handles,ax)
end



% --------------------------------------------------------------------
% POSITION CALCULATOR
% --------------------------------------------------------------------
function [rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles)

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Calculate the special positions
if strcmp(type,'columntype')
    if ~isequal(lenr,1)
        if b_isconstant(nrows)
            nrows = nrows(1);
            lenr = 1;
        else
            axno = sum(nrows);   % number of axes
            for i = 1:lenr
                cumm(i) = sum(nrows(1:i));
            end
            colend = cumm(1:end-1);
            colstart = colend + 1;
            rowend = [];
            rowstart = [];
        end
    end
    if lenr == 1
        axno = nrows * ncols;
        if ncols > 1
            colend = [1:ncols-1] * nrows;
            colstart = colend + 1;
        else
            colend = [];
            colstart = [];
        end
        if nrows > 1
            rowend = [(ncols-1)*nrows+1:nrows*ncols-1];
            rowstart = [2:nrows];
        else
            rowend = [];
            rowstart = [];
        end
    end
elseif strcmp(type,'rowtype')
    if ~isequal(lenc,1)
        if b_isconstant(ncols)
            ncols = ncols(1);
            lenc = 1;
        else
            axno = sum(ncols);   % number of axes
            for i = 1:lenc
                cumm(i) = sum(ncols(1:i));
            end
            rowend = cumm(1:end-1);
            rowstart = rowend + 1;
            colend = [];
            colstart = [];
        end
    end
    if lenc == 1
        axno = nrows * ncols;
        if ncols > 1
            colend = [(nrows-1)*ncols+1:nrows*ncols-1];
            colstart = [2:ncols];
        else
            colend = [];
            colstart = [];
        end
        if nrows > 1
            rowend = [1:nrows-1] * ncols;
            rowstart = rowend + 1;
        else
            rowend = [];
            rowstart = [];
        end
    end
end
handles.layout = [{nrows} {ncols}];
guidata(h,handles)



% --------------------------------------------------------------------
% TYPESWITCH
% --------------------------------------------------------------------
function typeswitch(h,handles)

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Get axes handles
A = handles.axes;

% Switch type
if strcmp(type,'rowtype')
    ncols = ncols(1);
    handles.type = 'columntype';
    NA = zeros(1,ncols*nrows);     % reorganise axes handles
    for i = 1:ncols
        for j = 1:nrows
            ind1 = (i - 1) * nrows + j;
            ind2 = (j - 1) * ncols + i;
            NA(ind1) = A(ind2);
        end
    end
    handles.axes = NA;
elseif strcmp(type,'columntype')
    nrows= nrows(1);
    handles.type = 'rowtype';
    NA = zeros(1,ncols*nrows);     % reorganise axes handles
    for i = 1:nrows
        for j = 1:ncols
            ind1 = (i - 1) * ncols + j;
            ind2 = (j - 1) * nrows + i;
            NA(ind1) = A(ind2);
        end
    end
    handles.axes = NA;
end
handles.layout = [{nrows} {ncols}];
guidata(h,handles)

% Refresh figure menus & context menus
for i = 1:length(handles.figure)
    refresh_figure_menus(h,handles,handles.figure(i))
    ax = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
    figure_cmenus(h,handles,ax)
end



% --------------------------------------------------------------------
% MENU REFRESHER
% --------------------------------------------------------------------
function menu_refresher(h,handles)

% Store default properties
default_handlevisibility = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','on');

% Calculate special positions
[rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles);
handles = guidata(h);

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
type = handles.type;
lenr = length(nrows);
lenc = length(ncols);

% Get axes handles
cax = gca;
A = handles.axes;
inx = find(A==cax);

% Define the menus
vmenu = handles.view_menu;
M = get(vmenu,'Children');
delete(M)

% Define callbacks for menu items
cb1 = 'b_ica_layouteditor';
cb2 = 'b_ica_gui4(''typeswitch'',gcbo,guidata(gcbo))';
cb3 = 'b_ica_gui4(''undock'',gcbo,guidata(gcbo))';

% Define the menu items
item1 = uimenu(vmenu,'Label','Layout','Callback',cb1);
if (isequal(lenc,1) & isequal(lenr,1)) | (lenr > 1 & b_isconstant(nrows)) | (lenc > 1 & b_isconstant(ncols))
    str = 'on';
else
    str = 'off';
end
item2 = uimenu(vmenu,'Label','Switch type','Callback',cb2,'Enable',str,'Accelerator','y');
if ~isempty(inx)
    item3 = uimenu(vmenu,'Label','Undock', 'Callback',cb3,'Separator','on','Accelerator','u');
end
if isequal(inx,1)
    if lenr == 1 & lenc == 1
        create_submenu_first(vmenu,'rowtype',nrows,ncols)      % menu item for first position
        create_submenu_first(vmenu,'columntype',nrows,ncols)      % menu item for first position
    else
        create_submenu_first(vmenu,type,nrows,ncols);      % menu item for first position
    end
end
if isequal(inx,length(A))
    if lenr == 1 & lenc == 1
        create_submenu_last(vmenu,'rowtype',nrows,ncols)      % menu item for last position
        create_submenu_last(vmenu,'columntype',nrows,ncols)      % menu item for last position
    else
        create_submenu_last(vmenu,type,nrows,ncols);      % menu item for last position
    end
end
if ~isempty(inx)
    if ~isempty(find(([0.5 rowend]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
        create_submenu_rowend(vmenu,nrows,ncols,inx)     % menu items for rowend position
    end
    if ~isempty(find(([0.5 rowstart]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
        create_submenu_rowstart(vmenu,nrows,ncols,inx)     % menu items for rowstart position
    end
    if ~isempty(find(([0.5 colend]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
        create_submenu_colend(vmenu,nrows,ncols,inx)     % menu items for colend position
    end
    if ~isempty(find(([0.5 colstart]-inx)==0))    % 0.5 is for avoiding empty == scalar comparison
        create_submenu_colstart(vmenu,nrows,ncols,inx)     % menu items for colstart position
    end
end

% Set children's uicontextmenus
for i = 1: length(A)
    ch = get(A(i),'Children');
    cm = get(A(i),'UIContextMenu');
    set(ch,'UIContextMenu',cm)
end

% Reset properties 
set(handles.figure1,'HandleVisibility',default_handlevisibility);

% --------------------------------------------------------------------
function varargout = refresh_figure_menus(h,handles,fig)

% Calculate the special positions
[rowstart,rowend,colstart,colend,axno] = position_calculator(h,handles);
handles = guidata(h);

% Get layout type
nrows = handles.layout{1};
ncols = handles.layout{2};
lenr = length(nrows);
lenc = length(ncols);

% Define the menus
vmenu = findobj(get(fig,'Children'),'Label','View');
delete(get(vmenu,'Children'))

% Define callbacks for 'View' menu items
posno = axno + 1;
cb = cell(1,posno+1);
cb{1} = 'closereq';
for i = 1:posno
    cb{i+1} = ['b_ica_gui4(''dock'',gcbo,guidata(gcbo),' num2str(i) ','''')'];
end

% Define the menu items for 'View' menu
flag{1} = 'st';
flag{2} = 'nd';
flag{3} = 'rd';
if posno > 3
    for i = 4:posno
        flag{i} = 'th';
    end
end

item = cell(1,posno+1);
item{1} = uimenu(vmenu,'Label','Close','Callback',cb{1},'Accelerator','q');
for i = 1:min(posno,9)    % only the 1-9 position has a shortcut
    item{i+1} = uimenu(vmenu,'Label',['Dock to ' num2str(i) flag{i} ' position'],...
        'Callback',cb{i+1},'Accelerator',[num2str(i)]);
end
if posno > 9 
    for i = 10:posno
        item{i+1} = uimenu(vmenu,'Label',['Dock to ' num2str(i) flag{i} ' position'],...
            'Callback',cb{i+1});
    end
end



% --------------------------------------------------------------------
% CREATE BUTTONS
% --------------------------------------------------------------------
function buts = create_buttons(h,handles,ax,button_struct)

% Check if the buttons allready exist
fig = gcf;
ch = get(fig,'Children');
button_handles = findobj(ch,'Style','pushbutton');
buttons = [];
for j = 1:length(button_handles)
    usd = get(button_handles(j),'UserData');
    if isequal(usd.axes_handle,ax)
        buttons(end+1) = button_handles(j);
    end
end
if isequal(length(buttons),length(button_struct))
    buts = [];
    return
else
    delete(buttons)
end

% Store default properties and set 'Units'
default_figure_units = get(fig,'Units');
default_axes_units = get(ax,'Units');
set(fig,'Units','pixels')
set(ax,'Units','pixels')

% Positions
ax_pos = get(ax,'Position');
lenb = length(button_struct);
buts = zeros(1,lenb);
for i = 1:lenb
    width = ax_pos(3) / lenb;
    left = ax_pos(1) + (i - 1) * width;
    buttom = ax_pos(2) + ax_pos(4);
    hight = 20;
    buttons_pos = [left buttom width hight];
    usd.axes_handle = ax;
    usd.button_number = lenb;
    buts(i) = uicontrol(fig,'Style','pushbutton','Position',buttons_pos,'String',...
        button_struct(i).string,'Callback',button_struct(i).callback,'UserData',usd);
end

% Reset properties 
set(fig,'Units',default_figure_units)
set(ax,'Units',default_axes_units)



% --------------------------------------------------------------------
% SET CURRENT AXES - callback for pushbuttons
% --------------------------------------------------------------------
function set_current_axes(h,handles,ax)

% Return if called from undocked figure
if isempty(handles)
    return
end

% Store default 'HandleVisibility'
default_handlevisibility = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','on');

% Set current axes
ch = get(handles.figure1,'Children');
A = findobj(ch,'Type','axes');
B = round(A*10000) / 10000;
inx = find(B==ax);
axes(A(inx))

% Set box color and refresh menus
Axes_ButtonDownFcn(h,[],handles)

% Reset 'HandleVisibility
set(handles.figure1,'HandleVisibility',default_handlevisibility);



% --------------------------------------------------------------------
% SET FONT WEIGHT of buttons - set current to 'bold'
% --------------------------------------------------------------------
function set_font_weight(h,handles,ax,varargin)

% Get current button handle
if ~isempty(varargin)
    curbut = varargin{1};
else
    curbut = gco;
end

% Get pushbutton handles
ch = get(gcf,'Children');
button_handles = findobj(ch,'Style','pushbutton');

% Set font weight
ax = round(ax*10000)/10000;
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    axes_handle = usd.axes_handle;
    axes_handle = round(axes_handle*10000)/10000;
    if isequal(axes_handle,ax)
        if isequal(button_handles(i),curbut)
            set(button_handles(i),'FontWeight','bold')
        else
            set(button_handles(i),'FontWeight','normal')
        end
    end
end



% --------------------------------------------------------------------
% GET CURRENT BUTTON
% --------------------------------------------------------------------
function fdi_index = get_current_button(h,handles,ax)

% Get pushbutton handles
ch = get(gcf,'Children');
button_handles = findobj(ch,'Style','pushbutton');
button_handles = flipud(button_handles);

% Get current substate ('fdi_index')
ax = round(ax*10000)/10000;
next = 1;
fdi_index = [];
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    axes_handle = usd.axes_handle;
    axes_handle = round(axes_handle*10000)/10000;
    if isequal(axes_handle,ax)
        fw = get(button_handles(i),'FontWeight');
        if strcmp(fw,'bold')
            fdi_index = next;
            break
        else
            next = next + 1;
        end
    end
end



% --------------------------------------------------------------------
% GET BUTTON HANDLE
% --------------------------------------------------------------------
function B = get_button_handle(h,handles,ax,no)

% Get pushbutton handles
ch = get(gcf,'Children');
button_handles = findobj(ch,'Style','pushbutton');
button_handles = flipud(button_handles);

% Get button handle of the 'no'-th button of 'ax'
ax = round(ax*10000)/10000;
ax_buttons = [];
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    axes_handle = usd.axes_handle;
    axes_handle = round(axes_handle*10000)/10000;
    if isequal(axes_handle,ax)
        ax_buttons(end+1) = button_handles(i);
    end
end

B = ax_buttons(no);



% --------------------------------------------------------------------
% DELETE BUTTONS
% --------------------------------------------------------------------
function delete_buttons(h,handles)

% Find buttons
button_handles = findobj(get(gcf,'Children'),'Style','pushbutton');

% Delete buttons
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    ax = usd.axes_handle;
    if isequal(ax,gca)
        delete(button_handles(i))
    end
end

% --------------------------------------------------------------------
function delete_buttons2(ax)

% Find buttons
button_handles = findobj(get(gcf,'Children'),'Style','pushbutton');

% Delete buttons
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    ax_handle = usd.axes_handle;
    if isequal(ax_handle,ax)
        delete(button_handles(i))
    end
end



% --------------------------------------------------------------------
% HANDLING THE LEGENDS
% --------------------------------------------------------------------
function [L,legend_pos] = legendcheck(h,ax)

% Initialize output
L = [];
legend_pos = [];

% Check if there is a legend for the axes whos handle is given in the input argument 'ax' 
ach = get(h,'Children');
legend_handles = findobj(ach,'Tag','legend');
if ~isempty(legend_handles)
    for i = 1:length(legend_handles)
        usd = get(legend_handles(i),'UserData');
        parent_axes = usd.PlotHandle;
        if isequal(ax,parent_axes)
            L = legend_handles(i);
            legend_pos = usd.legendpos;
            if length(legend_pos) > 1
                legend_pos = 1;         % default legend position
            end
            return
        end
    end
end

% --------------------------------------------------------------------
function lpos = legendposition(ha,legendpos,llen,lhgt)

%Get position vector from legendpos code
stickytol = 1;
cap = get(ha,'Position');
edge = 5;   % 5 point edge

if length(legendpos) == 4
    Pos = [legendpos(1) legendpos(2)+legendpos(4)-lhgt];
else
    switch legendpos
    case 0
        Pos = lscan(ha,llen,lhgt,0,stickytol);
    case 1
        Pos = [cap(1)+cap(3)-llen-edge cap(2)+cap(4)-lhgt-edge];
    case 2
        Pos = [cap(1)+edge cap(2)+cap(4)-lhgt-edge];
    case 3
        Pos = [cap(1)+edge cap(2)+edge];
    case 4
        Pos = [cap(1)+cap(3)-llen-edge cap(2)+edge];
    otherwise
        Pos = -1;
    end
end

if isequal(Pos,-1)
    lpos = [cap(1)+cap(3)-llen+edge cap(4)+cap(2)-lhgt llen lhgt];
else
    lpos = [Pos(1) Pos(2) llen lhgt];
end

% --------------------------------------------------------------------
function Pos = lscan(ha,wdt,hgt,tol,stickytol)

% Scan for good legend location.
% Calculate tile size
cap = get(ha,'Position');
xlim = get(ha,'Xlim');
ylim = get(ha,'Ylim');
H = ylim(2) - ylim(1);
W = xlim(2) - xlim(1);

dh = 0.03 * H;
dw = 0.03 * W;
Hgt = hgt * H / cap(4);
Wdt = wdt * W / cap(3);
Thgt = H / max(1,floor(H/(Hgt+dh)));
Twdt = W / max(1,floor(W/(Wdt+dw)));
dh = (Thgt - Hgt) / 2;
dw = (Twdt - Wdt) / 2;

% Get data, points and text
Kids = get(ha,'Children');
Xdata = [];
Ydata = [];
for i = 1:length(Kids),
    type = get(Kids(i),'Type');
    if strcmp(type,'line')
        xk = get(Kids(i),'XData');
        yk = get(Kids(i),'YData');
        n = length(xk);
        if n < 100 & n > 1
            xk = interp1(xk,linspace(1,n,200));
            yk = interp1(yk,linspace(1,n,200));
        end
        Xdata = [Xdata,xk];
        Ydata = [Ydata,yk];
    elseif strcmp(type,'patch') | strcmp(type,'surface')
        xk = get(Kids(i),'XData');
        yk = get(Kids(i),'YData');
        Xdata = [Xdata,xk(:)'];
        Ydata = [Ydata,yk(:)'];
    elseif strcmp(get(Kids(i),'type'),'text'),
        tmpunits = get(Kids(i),'Units');
        set(Kids(i),'Units','data')
        tmp = get(Kids(i),'Position');
        ext = get(Kids(i),'Extent');
        set(Kids(i),'Units',tmpunits);
        Xdata = [Xdata,[tmp(1) tmp(1)+ext(3)]];
        Ydata = [Ydata,[tmp(2) tmp(2)+ext(4)]];
    end
end
in = finite(Xdata) & finite(Ydata);
Xdata = Xdata(in);
Ydata = Ydata(in);

% Determine # of data points under each "tile"
xp = (0:Twdt:W-Twdt) + xlim(1);
yp = (0:Thgt:H-Thgt) + ylim(1);
wtol = Twdt / 100;
htol = Thgt / 100;
for j = 1:length(yp)
    if debug
        line([xlim(1) xlim(2)],[yp(j) yp(j)],'HandleVisibility','off');
    end
    for i = 1:length(xp)
        if debug
            line([xp(i) xp(i)],[ylim(1) ylim(2)],'HandleVisibility','off');
        end
        pop(j,i) = sum(sum((Xdata > xp(i)-wtol) & (Xdata < xp(i)+Twdt+wtol) & ...
            (Ydata > yp(j)-htol) & (Ydata < yp(j)+Thgt+htol)));    
    end
end

if all(pop(:) == 0)
    pop(1) = 1;
end

% Cover up fewest points.  After this while loop, pop will
% be lowest furthest away from the data
while any(pop(:) == 0)
    newpop = filter2(ones(3),pop);
    if all(newpop(:) ~= 0)
        break
    end
    pop = newpop;
end
[j,i] = find(pop == min(pop(:)));
xp =  xp - xlim(1) + dw;
yp =  yp - ylim(1) + dh;
Pos = [cap(1)+xp(i(end))*cap(3)/W
    cap(2)+yp(j(end))*cap(4)/H];

% --------------------------------------------------------------------
function delete_legend(h,handles)

% Find legend
[L,legendpos] = legendcheck(gcf,gca);

% Delete legend
if ~isempty(L)
    delete(L)
end

% --------------------------------------------------------------------
function delete_legend2(ax)

% Find legend
[L,legendpos] = legendcheck(gcf,ax);

% Delete legend
if ~isempty(L)
    delete(L)
end



% --------------------------------------------------------------------
% UNDOCKED FIGURE RESIZE FUNCTION
% --------------------------------------------------------------------
function resize_undocked

% Get button handles
ch = get(gcf,'Children');
button_handles = findobj(ch,'Style','pushbutton');
button_handles = flipud(button_handles);

% Get axes handles
A = findobj(ch,'Type','axes','Tag','');

% Move buttons
if ~isempty(button_handles)     % remember old 'Units'
    old_units_buttons = get(button_handles(1),'Units');
end
old_axes_units = get(A(1),'Units');
set(button_handles,'Units','pixels');   % set 'Units'
set(A,'Units','pixels');
fig_pos = get(gcf,'Position');

next = 1;
ax = 0;
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    if ~isequal(usd.axes_handle,ax)
        next = 1;
    end
    ax = usd.axes_handle;
    default_units = get(ax,'Units');
    set(ax,'Units','pixels')
    lenb = usd.button_number;
    ax_pos = get(ax,'Position');
    width = ax_pos(3) / lenb;
    left = ax_pos(1) + (next - 1) * width;
    buttom = ax_pos(2) + ax_pos(4);
    hight = 20;
    if ax_pos(2) + ax_pos(4) + hight > fig_pos(4)     % decrease hight if button exceeds figure window
        ax_hight = fig_pos(4) - (ax_pos(2) + hight);
        set(ax,'Position',[ax_pos(1) ax_pos(2) ax_pos(3) ax_hight])
        buttom = ax_pos(2) + ax_hight;
    end
    buttons_pos = [left buttom width hight];
    set(button_handles(i),'Position',buttons_pos);
    set(ax,'Units',default_units)
    next = next + 1;
end

if ~isempty(button_handles)
    set(button_handles,'Units',old_units_buttons);   % reset 'Units'
end
set(A,'Units',old_axes_units)



% --------------------------------------------------------------------
% UNDOCKED FIGURE CLOSE REQUEST FUNCTION
% --------------------------------------------------------------------
function closereq_undocked

try
    cru
catch
    delete(gcf)
end

% --------------------------------------------------------------------
function cru

% Get 'handles' structure
global ICA_GUI_HANDLE
h = ICA_GUI_HANDLE;
handles = guidata(h);

% Delete figure handle from 'handles' structure
cax = gca;
for i = 1:length(handles.figure)
    A(i) = findobj(get(handles.figure(i),'Children'),'Type','axes');
end
inx = find(A==cax);
handles.figure(inx) = [];
handles.figdisplay(inx) = [];
guidata(h,handles)

% Delete figure
delete(gcf)



% --------------------------------------------------------------------
% UNDOCKED FIGURE'S ZOOM SUBMENU CALLBACK
% --------------------------------------------------------------------
function figure_zoom

% Set 'Checked' property and switch zoom state
zoom_submenu = findobj(get(gcf,'Children'),'Label','Zoom');
isch = get(zoom_submenu,'Checked');
if strcmp(isch,'on')
    set(zoom_submenu,'Checked','off')
    zoom off
else
    set(zoom_submenu,'Checked','on')
    zoom on
end



% --------------------------------------------------------------------
% SET UNDOCKED FIGURE'S KEYPRESS FUNCTION
% --------------------------------------------------------------------
function set_keypressfcn(h,handles)

% Get axes handles
cax = gca;
A = handles.axes;
B = 0;
for i = 1:length(handles.figure)
    B(i) = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
end
inx = find(B==cax);

% Set keypress function
if ~isempty(inx)
    kpf = handles.figdisplay(inx).keypress;
    set(gcf,'KeypressFcn',kpf)
end



% --------------------------------------------------------------------
% SET UNDOCKED FIGURE'S ZOOM
% --------------------------------------------------------------------
function set_zoom(h,handles)

% Get handle
ch = get(gcf,'Children');
zo = findobj(ch,'Label','Zoom');
isch = get(zo,'Checked');

% Set zoom state
if strcmp(isch,'on')
    zoom on
else
    zoom off
end



% --------------------------------------------------------------------
% SET 'DISPLAY' FIELD
% --------------------------------------------------------------------
function set_display_field(h,handles,str,varargin)

% Get axes handles and decide whether current axes is on GUI or on undocked figure
A = handles.axes;
cax = gca;
inx = find(A==cax);
if ~isempty(inx)
    figorgui = 'gui';
else
    clear A
    for i = 1:length(handles.figure)
        A(i) = findobj(get(handles.figure(i),'Children'),'Type','axes','Tag','');
    end
    inx = find(A==cax);
    figorgui = 'fig';
end

% Set 'display' field
display.name = str;
if strcmp(str,'variance_plot')      % 'settings' field
    display.settings = {getappdata(cax,'clustering_method'),...
            get(handles.inter_submenu,'Checked'),get(handles.intra_submenu,'Checked'),...
            get(handles.extra_submenu,'Checked'),get(handles.first_submenu,'Checked'),...
            get(handles.allfirst_submenu,'Checked'),get(handles.hold_submenu,'Checked'),...
            get(handles.legend1_submenu,'Checked')};
elseif strcmp(str,'structural_elements')
    display.settings = {getappdata(cax,'clustering_method'),...
            get(handles.superimpose1_submenu,'Checked'),get(handles.superimpose2_submenu,'Checked'),...
            get(handles.superimpose3_submenu,'Checked'),get(handles.superimpose4_submenu,'Checked'),...
            get(handles.crossline_submenu,'Checked')};
elseif strcmp(str,'wavelet_average')
    display.settings = {get(handles.hz1_3_submenu,'Checked'),...
            get(handles.hz3_6_submenu,'Checked'),get(handles.hz6_20_submenu,'Checked'),...
            get(handles.hz20_50_submenu,'Checked'),get(handles.legend2_submenu,'Checked')};
elseif ~isempty(find(strcmp(str,{'return_plot','burst_parameters','lomb_periodogram',...
            'burstwave_matrices','scale_cv_functions'})))
    display.settings = {getappdata(cax,'clustering_method')};
else
    display.settings = {};
end
curbut = get_current_button(h,handles,cax);        % 'button' field
display.button = curbut;
if ~isempty(find(strcmp(str,{'wavelet','burstwave_matrices','scale_cv_functions',...        % 'stringinput' field
            'spike_triggered_wavelet','spike_triggered_wavelet_average',...
            'spstw','spstw_average'})))
    display.stringinput = varargin{1};
else
    display.stringinput = '';
end
if strcmp(str,'structural_elements')        % 'keypress' field
    display.keypress = 'b_ica_gui4(''callgreyplot'')';
else
    display.keypress = '';
end
if strcmp(str,'burst_parameters')       % 'opening button' field
    display.opening_button = 'current';
else
    display.opening_button = 1;
end

if strcmp(figorgui,'fig')
    handles.figdisplay(inx) = display;
elseif strcmp(figorgui,'gui')
    handles.guidisplay(inx) = display;
end
guidata(h,handles)



% --------------------------------------------------------------------
% SET AXES PROPERTIES
% --------------------------------------------------------------------
function set_axes_properties(h,handles)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
end

% Set axes properties
ch = get(gca,'Children');
uicm = get(gca,'UIContextMenu');
set(ch,'UIContextMenu',uicm)
if isequal(gcf,handles.figure1)
    bdf = 'b_ica_gui4(''Axes_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes,'XColor','black','YColor','black')
    set(gca,'XColor','red','YColor','red','ButtonDownFcn',bdf)
    set(ch,'ButtonDownFcn',bdf)
end



% --------------------------------------------------------------------
% AXES BUTTONDOWN FUNCTION
% --------------------------------------------------------------------
function varargout = Axes_ButtonDownFcn(h,eventdata,handles,varargin)

% Set box color
set(handles.axes,'XColor','black')
set(handles.axes,'YColor','black')
set(gca,'XColor','red')
set(gca,'YColor','red')

% Refresh 'View' menu
menu_refresher(h,handles);



% --------------------------------------------------------------------
% CLOSE REQUEST FUNCTION
% --------------------------------------------------------------------
function varargout = figure1_CloseRequestFcn(h, eventdata, handles, varargin)

try
% Delete figures
    delete(handles.figure(find(ishandle(handles.figure))))

% Delete GUI
    delete(handles.figure1)
    
catch
    delete(gcf)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             O P E N                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------------------------------------------------------------
% CELL SELECTION
% --------------------------------------------------------------------
function varargout = open_submenu_Callback(h,eventdata,handles,varargin)

% Launch Open GUI
b_ica_open



% --------------------------------------------------------------------
% NEXT CELL
% --------------------------------------------------------------------
function varargout = next_submenu_Callback(h,eventdata,handles,varargin)

% Get open GUI handle
global OPEN_GUI_HANDLE
g = OPEN_GUI_HANDLE;
handles2 = guidata(g);

% Set listbox 'Value'
index_selected = get(handles2.listbox1,'Value');
file_list = get(handles2.listbox1,'String');
if length(file_list) < index_selected+1
    return
end
set(handles2.listbox1,'Value',index_selected+1)

% Open
b_ica_open('listbox1_Callback',g,[],handles2);



% --------------------------------------------------------------------
% PREVIOUS CELL
% --------------------------------------------------------------------
function varargout = previous_submenu_Callback(h,eventdata,handles,varargin)

% Get open GUI handle
global OPEN_GUI_HANDLE
g = OPEN_GUI_HANDLE;
handles2 = guidata(g);

% Set listbox 'Value'
index_selected = get(handles2.listbox1,'Value');
if index_selected-1 < 1
    return
end
set(handles2.listbox1,'Value',index_selected-1)

% Open
b_ica_open('listbox1_Callback',g,[],handles2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             F I L L                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------------------------------------------------------------
% CONTEXT MENU CALLBACKS
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% VARIANCE PLOT
% --------------------------------------------------------------------
function variance_plot(h,handles,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

% Get Clustering menu settings if called from open GUI
if ~isempty(varargin)
    method = varargin{1};
    isinter = varargin{2};
    isintra = varargin{3};
    isextra = varargin{4};
    isfirst = varargin{5};
    isallfirst = varargin{6};
    isholdchecked = varargin{7};
    islegend = varargin{8};
else

% Get Clustering menu settings if called from ICA GUI
    method = handles.clustering_method;
    isintra = get(handles.intra_submenu,'Checked');
    isextra = get(handles.extra_submenu,'Checked');
    isinter = get(handles.inter_submenu,'Checked');
    isfirst = get(handles.first_submenu,'Checked');
    isallfirst = get(handles.allfirst_submenu,'Checked');
    isholdchecked = get(handles.hold_submenu,'Checked');
    islegend = get(handles.legend1_submenu,'Checked');
end

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Load data
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)

% Plot
d = length(IntraBurstIvCv);
plot_handle = [];
legend_string = {};
plot(1,'Color','white')         % this line garantates that current axes is not the legend
delete(get(gca,'Children'))
if strcmp(isintra,'on')
    plot_handle(end+1) = plot([1:d],IntraBurstIvCv,'g');
    hold on
    legend_string{end+1} = 'intraburstivcv';
end
if strcmp(isextra,'on')
    plot_handle(end+1) = plot([1:d],ExtraBurstIvCv,'b');
    hold on
    legend_string{end+1} = 'extraburstivcv';
end
if strcmp(isinter,'on')
    plot_handle(end+1) = plot([1:d],InterBurstIvCv,'c');
    hold on
    legend_string{end+1} = 'interburstivcv';
end
if strcmp(isfirst,'on')
    plot_handle(end+1) = plot([1:d],FirstSpikeCv,'m');
    hold on
    legend_string{end+1} = 'firstspikecv';
end
if strcmp(isallfirst,'on')
    plot_handle(end+1) = plot([1:d],AllFirstSpikeCv,'k');
    hold on
    legend_string{end+1} = 'allfirstspikecv';
end
if strcmp(islegend,'on')
    if ~isempty(legend_string)
        L = legend(legend_string,1);
        set(L,'Units','characters');
    end
end
xlim([1 d]);
x_lim = xlim;
y_lim = ylim;
kk = (y_lim(2) - y_lim(1)) / 3;
length_vdisc = length(Vdisc);
if strcmp(isholdchecked,'off')
    if length_vdisc < 100,
        T0 = text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length_vdisc)],'Color','r');
    else
        T0 = text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length_vdisc)],'Color','k');
    end
    set(T0,'Units','normalized');
    hold off
end

T1 = text(x_lim(1)+(x_lim(2)-x_lim(1))/2,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,...
    ['{\itCcc : }',num2str(Ccc)],'Color','k');      % Cophenetic coefficient
T2 = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,...
    [method],'Color','k');      % Method
set(T1,'Units','normalized');
set(T2,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'variance_plot')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end

% Refresh 'View' menu
if ~isundocked
    menu_refresher(h,handles);
end



% --------------------------------------------------------------------
% STRUCTURAL ELEMENTS
% --------------------------------------------------------------------
function structural_elements(h,handles,fdi_index,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Create 'time' vector
dt = 0.0001;
len = i_second - i_first;
time = [1:i_second] * dt; 

% Find the different states
fdi = find(diff(IntraBurstIvCv));
fdi = fdi + 1;
lenfdi = length(fdi);
dec = fdi(fdi_index);

% Create buttons
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''structural_elements'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''structural_elements'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Plot
z = zeros(1,i_second);
z(Vdisc) = 1;
z(Vdisc+1) = -1;
hold off
plot_handle = plot(time(i_first:i_second),z(i_first:i_second),'b');
axis([time(i_first) time(i_second) -1.5 1.5])
hold on;
sbd = size(Burst{dec},2);
hans = zeros(1,sbd);
for j = 1:sbd
    ind1 = Vdisc(Burst{dec}(1,j));
    ind2 = Vdisc(Burst{dec}(2,j));
    hans(j) = plot(time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r');
end
hold off

if ~isempty(varargin)       % get settings if called from open GUI
    is1 = varargin{2};
    is2 = varargin{3};
    is3 = varargin{4};
    is4 = varargin{5};
    isc = varargin{6};
else                        % get settings if called from ICA GUI
    is1 = get(handles.superimpose1_submenu,'Checked');
    is2 = get(handles.superimpose2_submenu,'Checked');
    is3 = get(handles.superimpose3_submenu,'Checked');
    is4 = get(handles.superimpose4_submenu,'Checked');
    isc = get(handles.crossline_submenu,'Checked');
end

if strcmp(is1,'on')
    superimpose1(h,handles);
end
if strcmp(is2,'on')
    superimpose2(h,handles);
end
if strcmp(is3,'on')
    superimpose3(h,handles);
end
if strcmp(is4,'on')
    superimpose4(h,handles);
end
if strcmp(isc,'on')
    crossline(h,handles);
end

x_lim = xlim;
y_lim = ylim;
T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'line_handles' field and 'greyplot' fields to empty matrix
handles.line_handles = [];
handles.greyplot1_handle = [];
handles.greyplot2_handle = [];
handles.greyplot3_handle = [];
handles.greyplot4_handle = [];
guidata(h,handles)

% Set 'display' field
set_display_field(h,handles,'structural_elements')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end

% --------------------------------------------------------------------
function superimpose1(h,handles)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
end

% Load wavelet average vectors
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

sw = size(Wavevec,2);
ttm = linspace(i_first/10000,i_second/10000,sw);

% Plot wavelet average magnitudes from 1 to 3 hz against time on the spike train
hold on
if isempty(handles.greyplot1_handle)
    wm = Wavevec(4,:) - mean(Wavevec(4,:));
    wwm = wm / max(wm) * 1.5;
    handles.greyplot1_handle = plot(ttm,wwm,'k');
    if ~isempty(handles.greyplot2_handle)
        set(handles.greyplot2_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot3_handle)
        set(handles.greyplot3_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot4_handle)
        set(handles.greyplot4_handle,'Visible','off')
    end
else
    vis = get(handles.greyplot1_handle,'Visible');
    if strcmp(vis,'on')
        set(handles.greyplot1_handle,'Visible','off')
    else
        set(handles.greyplot1_handle,'Visible','on')
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    end
end
hold off
guidata(h,handles)

% Set axes properties
set_axes_properties(h,handles)

% --------------------------------------------------------------------
function superimpose2(h,handles)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
end

% Load wavelet average vectors
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

sw = size(Wavevec,2);
ttm = linspace(i_first/10000,i_second/10000,sw);

% Plot wavelet average magnitudes from 3 to 6 hz against time on the spike train
hold on
if isempty(handles.greyplot2_handle)
    wm = Wavevec(3,:) - mean(Wavevec(3,:));
    wwm = wm / max(wm) * 1.5;
    handles.greyplot2_handle = plot(ttm,wwm,'k');
    if ~isempty(handles.greyplot1_handle)
        set(handles.greyplot1_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot3_handle)
        set(handles.greyplot3_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot4_handle)
        set(handles.greyplot4_handle,'Visible','off')
    end
else
    vis = get(handles.greyplot2_handle,'Visible');
    if strcmp(vis,'on')
        set(handles.greyplot2_handle,'Visible','off')
    else
        set(handles.greyplot2_handle,'Visible','on')
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    end
end
hold off
guidata(h,handles)

% Set axes properties
set_axes_properties(h,handles)

% --------------------------------------------------------------------
function superimpose3(h,handles)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
end

% Load wavelet average vectors
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

sw = size(Wavevec,2);
ttm = linspace(i_first/10000,i_second/10000,sw);

% Plot wavelet average magnitudes from 6 to 20 hz against time on the spike train
hold on
if isempty(handles.greyplot3_handle)
    wm = Wavevec(2,:) - mean(Wavevec(2,:));
    wwm = wm / max(wm) * 1.5;
    handles.greyplot3_handle = plot(ttm,wwm,'k');
    if ~isempty(handles.greyplot1_handle)
        set(handles.greyplot1_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot2_handle)
        set(handles.greyplot2_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot4_handle)
        set(handles.greyplot4_handle,'Visible','off')
    end
else
    vis = get(handles.greyplot3_handle,'Visible');
    if strcmp(vis,'on')
        set(handles.greyplot3_handle,'Visible','off')
    else
        set(handles.greyplot3_handle,'Visible','on')
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    end
end
hold off
guidata(h,handles)

% Set axes properties
set_axes_properties(h,handles)

% --------------------------------------------------------------------
function superimpose4(h,handles)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
end

% Load wavelet average vectors
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

sw = size(Wavevec,2);
ttm = linspace(i_first/10000,i_second/10000,sw);

% Plot wavelet average magnitudes from 20 to 50 hz against time on the spike train
hold on
if isempty(handles.greyplot4_handle)
    wm = Wavevec(1,:) - mean(Wavevec(1,:));
    wwm = wm / max(wm) * 1.5;
    handles.greyplot4_handle = plot(ttm,wwm,'k');
    if ~isempty(handles.greyplot1_handle)
        set(handles.greyplot1_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot2_handle)
        set(handles.greyplot2_handle,'Visible','off')
    end
    if ~isempty(handles.greyplot3_handle)
        set(handles.greyplot3_handle,'Visible','off')
    end
else
    vis = get(handles.greyplot4_handle,'Visible');
    if strcmp(vis,'on')
        set(handles.greyplot4_handle,'Visible','off')
    else
        set(handles.greyplot4_handle,'Visible','on')
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
    end
end
hold off
guidata(h,handles)

% Set axes properties
set_axes_properties(h,handles)

% --------------------------------------------------------------------
function crossline(h,handles)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
end
method = getappdata(gca,'clustering_method');

% Load data
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Create 'time' vector
dt = 0.0001;
len = i_second - i_first;
time = [1:i_second] * dt; 

% Get pushbutton handles
ch = get(gcf,'Children');
button_handles = findobj(ch,'Style','pushbutton');
button_handles = flipud(button_handles);

% Get current substate ('fdi_index')
ax = gca;
fdi_index = get_current_button(h,handles,ax);

% Find the different states
fdi = find(diff(IntraBurstIvCv));
fdi = fdi + 1;
lenfdi = length(fdi);
dec = fdi(fdi_index);

% Drawing the lines or switching their visibility
hold on;
sbd = size(Burst{dec},2);
if isempty(handles.line_handles)
    line_handles = zeros(1,sbd);
    for j = 1:sbd
        rajz = [time(Vdisc(Burst{dec}(1,j))) time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
        line_handles(j) = line(rajz(1,:),rajz(2,:),'Color','k');
    end
    handles.line_handles = line_handles;
else
    vis = get(handles.line_handles,'Visible');
    if strcmp(vis{1},'on')
        set(handles.line_handles,'Visible','off')
    else
        set(handles.line_handles,'Visible','on')
    end
end
hold off
guidata(h,handles)

% Set axes properties
set_axes_properties(h,handles)



% --------------------------------------------------------------------
% RETURN PLOTS
% --------------------------------------------------------------------
function return_plots(h,handles,fdi_index,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Find the different states
fdi = find(diff(IntraBurstIvCv));
fdi = fdi + 1;
lenfdi = length(fdi);
dec = fdi(fdi_index);

% Create buttons
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''return_plots'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''return_plots'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Plot
ivs = diff(Vdisc);
sbd = size(Burst{dec},2);
intra_xdata = [];
intra_ydata = [];
for j = 1:sbd
    ind1 = Burst{dec}(1,j) + 1;
    ind2 = Burst{dec}(2,j) - 1;
    intra_xdata = [intra_xdata ivs(ind1-1:ind2-1)];
    intra_ydata = [intra_ydata ivs(ind1:ind2)];
end

hold off
plot_handle1 = plot(ReturnPlotXData,ReturnPlotYData,'.');
hold on
plot_handle2 = plot(intra_xdata,intra_ydata,'.','Color','green');
hold off

x_lim = xlim;
y_lim = ylim;
T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+4*(y_lim(2)-y_lim(1))/5,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'return_plots')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% BURST PARAMETERS
% --------------------------------------------------------------------
function burst_parameters(h,handles,param,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_PARAM_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Create buttons
input_string{1} = 'Burstiness';
input_string{2} = 'IntraBurstFrequency';
input_string{3} = 'IntraBurstSpikeNumber';
input_string{4} = 'BurstLength';
input_string{5} = 'BurstFrequency';
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:5
        spec = ['b_ica_gui4(''burst_parameters'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = input_string{i};
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:5
        spec = ['b_ica_gui4(''burst_parameters'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = input_string{i};
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Find the different states
fdi = find(diff(IntraBurstFrequency));
fdi = fdi + 1;
lenfdi = length(fdi);

% Plot
hold off
switch param
case 1
    b_bar2('blue',Burstiness(fdi))
    x_lim = xlim;
    y_lim = ylim;
    T1 = text(x_lim(1)+3*(x_lim(2)-x_lim(1))/4,y_lim(1)+4*(y_lim(2)-y_lim(1))/5,'Burstiness',...
        'Color','blue','HorizontalAlignment','center');
case 2
    b_bar2('red',IntraBurstFrequency(fdi))
    x_lim = xlim;
    y_lim = ylim;
    T1 = text(x_lim(1)+3*(x_lim(2)-x_lim(1))/4,y_lim(1)+4*(y_lim(2)-y_lim(1))/5,'IntraBurstFrequency',...
        'Color','red','HorizontalAlignment','center');
case 3
    b_bar2('green',IntraBurstSpikeNumber(fdi))
    x_lim = xlim;
    y_lim = ylim;
    T1 = text(x_lim(1)+3*(x_lim(2)-x_lim(1))/4,y_lim(1)+4*(y_lim(2)-y_lim(1))/5,'IntraBurstSpikeNumber',...
        'Color','green','HorizontalAlignment','center');
case 4
    b_bar2('cyan',BurstLength(fdi))
    x_lim = xlim;
    y_lim = ylim;
    T1 = text(x_lim(1)+3*(x_lim(2)-x_lim(1))/4,y_lim(1)+4*(y_lim(2)-y_lim(1))/5,'BurstLength',...
        'Color','cyan','HorizontalAlignment','center');
case 5
    b_bar2('yellow',BurstFrequency(fdi))
    x_lim = xlim;
    y_lim = ylim;
    T1 = text(x_lim(1)+3*(x_lim(2)-x_lim(1))/4,y_lim(1)+4*(y_lim(2)-y_lim(1))/5,'BurstFrequency',...
        'Color','yellow','HorizontalAlignment','center');
end

x_lim = xlim;
y_lim = ylim;
T2 = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T1,'Units','normalized');
set(T2,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'burst_parameters')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% INSTANTENOUS FREQUENCY
% --------------------------------------------------------------------
function instantenous_frequency(h,handles,long_index,varargin)

% Define long_index if called from OPEN GUI
if exist('long_index') == 0
    long_index = 1;
end

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\thres'];
dd = dir(pth);
kk = zeros(1,length(dd));
for i = 3:length(dd)
    if ~dd(i).isdir
        kk(i) = strcmp(dd(i).name(17:22),CELL(1:6));
    else
        kk(i) = logical(0);
    end
end
nn = find(kk);
fln = dd(nn(long_index)).name;
needbuts = 0;
if length(nn) == 2
    needbuts = 1;
elseif length(nn) > 2
    error('Long registration.')
end
ff = [pth '\' fln];
load(ff)
fs = findstr(fln,'_');
fname = fln(fs(3)+1:fs(5)-1);
datinx1 = str2num(fln(fs(5)+1:fs(6)-1));
datinx2 = str2num(fln(fs(6)+1:end-4));
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Create buttons
if needbuts
    sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    if isequal(gcf,handles.figure1)
        for i = 1:2
            spec = ['b_ica_gui4(''spstw'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw sca spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    else
        for i = 1:2
            spec = ['b_ica_gui4(''spstw'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    end
    if ~isempty(buts)
        fw = findobj(buts,'FontWeight','bold');
        if isempty(fw)
            set_font_weight(h,handles,gca,buts(1))
        end
    end
end

% Calculate instantenous frequency
dt = 0.0001;
len = datinx2 - datinx1;
isi = diff(vdisc) * dt;
instfrek = zeros(1,len);
for i = 1:length(vdisc)-1
    instfrek(vdisc(i):vdisc(i+1)) = 1 / isi(i);
end
instfrek(1:vdisc(1)-1) = 1 / (vdisc(1) * dt);
instfrek(vdisc(end):len) = 1 / ((len - vdisc(end)) * dt);
    
% Create 'time' vector
dt = 0.0001;
len = datinx2 - datinx1;
time = [1:len] * dt;

% Plot
instfrek2 = instfrek(i_first:i_second);
time2 = time(i_first:i_second);
hold off
plot(time2,instfrek2);
y_lim = ylim;
axis([time2(1) time2(end) y_lim(1) y_lim(2)])

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'instantenous_frequency')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% LOMB PERIODOGRAM
% --------------------------------------------------------------------
function lomb_periodogram(h,handles,fdi_index,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_LOMB_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Create buttons
lenfdi = length(Pxx);
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''lomb_periodogram'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''lomb_periodogram'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Plot Lomb-spectrum and significance levels
z05 = Z{fdi_index}(1);
z02 = Z{fdi_index}(2);
z01 = Z{fdi_index}(3);
z001 = Z{fdi_index}(4);

hold off
period_handle = fill(Pxx{fdi_index},Pyy{fdi_index},'r');
y_lim = ylim;
axis([0 max(Pxx{fdi_index}) y_lim(1) y_lim(2)]) % x: frequency [Hz]; y: PSD [s^2]
hold on
v = axis;
plot([v(1) v(2)], [ z05 z05], 'k:')     % 5% significance
plot([v(1) v(2)], [ z02 z02], 'k:')     % 2% significance
plot([v(1) v(2)], [ z01 z01], 'k:')     % 1% significance
plot([v(1) v(2)], [ z001 z001], 'k:')   % 0.1% significance
hold off

x_lim = xlim;
y_lim = ylim;
T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'lomb_periodogram')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% WAVELET
% --------------------------------------------------------------------
function wavelet(h,handles,strinp,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

if strcmp(strinp,'power')
    st = 'WAVELET_POWER_';
elseif strcmp(strinp,'phase')
    st = 'WAVELET_PHASE_';
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavelet\' strinp];
fln = [st CELL '.jpg'];
ff = [pth '\' fln];
str = ['WaveletImage = imread(''' ff ''');'];
eval(str)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Plot
hold off
imagesc(WaveletImage)

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'wavelet',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% WAVELET AVERAGE
% --------------------------------------------------------------------
function wavelet_average(h,handles,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Get settings if called from open GUI
if ~isempty(varargin)
    is1_3 = varargin{1};
    is3_6 = varargin{2};
    is6_20 = varargin{3};
    is20_50 = varargin{4};
    islegend = varargin{5};
else

% Get settings if called from ICA GUI
    is1_3 = get(handles.hz1_3_submenu,'Checked');
    is3_6 = get(handles.hz3_6_submenu,'Checked');
    is6_20 = get(handles.hz6_20_submenu,'Checked');
    is20_50 = get(handles.hz20_50_submenu,'Checked');
    islegend = get(handles.legend2_submenu,'Checked');
end

% Plot
sw = size(Wavevec,2);
ttm = linspace(i_first/10000,i_second/10000,sw);

legend_string = {};
delete(get(gca,'Children'))
hold off
if strcmp(is20_50,'on')
    plot1_handle = plot(ttm,Wavevec(1,:),'r');
    hold on
    legend_string{end+1} = '20 - 50 Hz';
end
if strcmp(is6_20,'on')
    plot1_handle = plot(ttm,Wavevec(2,:),'b');
    hold on
    legend_string{end+1} = '6 - 20 Hz';
end
if strcmp(is3_6,'on')
    plot1_handle = plot(ttm,Wavevec(3,:),'g');
    hold on
    legend_string{end+1} = '3 - 6';
end
if strcmp(is1_3,'on')
    plot1_handle = plot(ttm,Wavevec(4,:),'m');
    hold on
    legend_string{end+1} = '1 - 3 Hz';
end
if strcmp(islegend,'on')
    if ~isempty(legend_string)
        L = legend(legend_string,2);
        set(L,'Units','characters');
    end
end
hold off
y_lim = ylim;
intlen =  length(Wavevec);
axis([ttm(1) ttm(end) y_lim(1) y_lim(2)]);

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'wavelet_average')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% WAVELET AVERAGE FFT
% --------------------------------------------------------------------
function wavelet_average_fft(h,handles,no)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

ff = [DATAPATH 'ICA\ica_newave3\wavelet_settings\newstep.mat'];
load(ff)
samprate = 10000 / Newstep;

% Create buttons
input_string{1} = '1 - 3 Hz';
input_string{2} = '3 - 6 Hz';
input_string{3} = '6 - 20 Hz';
input_string{4} = '20 - 50 Hz';
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:4
        spec = ['b_ica_gui4(''wavelet_average_fft'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = input_string{i};
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:4
        spec = ['b_ica_gui4(''wavelet_average_fft'',gcbo,guidata(gcbo),' num2str(i) ')'];
        button_struct(i).string = input_string{i};
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Plot
hold off
lenwv = size(Wavevec,2);
switch no
case 1
    [p_1 f_1] = b_fft2(Wavevec(4,:),samprate,lenwv,2);
    p_1 = p_1';
    f_1 = f_1';
    fnd = find(f_1<5);
    lim1 = fnd(end);
    p_11 = [0;p_1(2:lim1);0];
    f_11 = [0;f_1(2:lim1);f_1(lim1)];
    fill(f_11,p_11,'m');
    ym = max([p_1(2:lim1)]);
    x_lim = xlim;
    xind = x_lim(1) + ((x_lim(2) - x_lim(1)) / 2); 
    T = text(xind,0.5*ym,'{\it1 - 3 Hz}');
case 2
    [p_2 f_2] = b_fft2(Wavevec(3,:),samprate,lenwv,2);
    p_2 = p_2';
    f_2 = f_2';
    fnd = find(f_2<5);
    lim2 = fnd(end);
    p_22 = [0;p_2(2:lim2);0];
    f_22 = [0;f_2(2:lim2);f_2(lim2)];
    fill(f_22,p_22,'g');
    ym = max([p_2(2:lim2)]);
    x_lim = xlim;
    xind = x_lim(1) + ((x_lim(2) - x_lim(1)) / 2); 
    T = text(xind,0.5*ym,'{\it3 - 6 Hz}');
case 3
    [p_3 f_3] = b_fft2(Wavevec(2,:),samprate,lenwv,2);
    p_3 = p_3';
    f_3 = f_3';
    fnd = find(f_3<10);
    lim3 = fnd(end);
    p_33 = [0;p_3(2:lim3);0];
    f_33 = [0;f_3(2:lim3);f_3(lim3)];
    fill(f_33,p_33,'b');
    ym = max([p_3(2:lim3)]);
    x_lim = xlim;
    xind = x_lim(1) + ((x_lim(2) - x_lim(1)) / 2); 
    T = text(xind,0.5*ym,'{\it6 - 20 Hz}');
case 4
    [p_4 f_4] = b_fft2(Wavevec(1,:),samprate,lenwv,2);
    p_4 = p_4';
    f_4 = f_4';
    fnd = find(f_4<20);
    lim4 = fnd(end);
    p_44 = [0;p_4(2:lim4);0];
    f_44 = [0;f_4(2:lim4);f_4(lim4)];
    fill(f_44,p_44,'r');
    ym = max([p_4(2:lim4)]);
    x_lim = xlim;
    xind = x_lim(1) + ((x_lim(2) - x_lim(1)) / 2); 
    T = text(xind,0.5*ym,'{\it20 - 50 Hz}');
end
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'wavelet_average_fft')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% BURSTWAVE MATRICES
% --------------------------------------------------------------------
function burstwave_matrices(h,handles,fdi_index,strinp,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method '\' strinp];
fln = ['BURSTWAVE_' method '_' CELL '_substate' num2str(fdi_index) '.jpg'];
ff = [pth '\' fln];
str = ['BurstWaveImage = imread(''' ff ''');'];
eval(str)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Find the different states
fdi = find(diff(IntraBurstIvCv));
fdi = fdi + 1;
lenfdi = length(fdi);

% Create buttons
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''burstwave_matrices'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''burstwave_matrices'',gcbo,guidata(gcbo),' strinp ',' num2str(i) ',''' strinp ''')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Plot
hold off
imagesc(BurstWaveImage)
y_lim = ylim;
x_lim = xlim;

T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'burstwave_matrices',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% SCALE-CV FUNCTIONS
% --------------------------------------------------------------------
function scale_cv_functions(h,handles,strinp,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_EEG_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));
if strcmp(strinp,'power')
    ScaleCv = ScaleCvPower;
elseif strcmp(strinp,'phase')
    ScaleCv = ScaleCvPhase;
end

global DATAPATH
ff = [DATAPATH 'ICA\ica_newave3\wavelet_settings\f.mat'];
load(ff)

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Plot
lenfdi = length(ScaleCv);
hold off
clr = zeros(1,lenfdi);
for b = 1:lenfdi
    clr(b) = (b - 1) / (lenfdi - 1);
    semilogx(ScaleVector,ScaleCv{b},'Color',[0 1-clr(b) clr(b)]);
    hold on
    x_lim = xlim;
    y_lim = ylim;
    axis([1 ScaleVector(1) y_lim(1) y_lim(2)]);
    legend_matrix{b} = ['substate: ' num2str(b) '/' num2str(lenfdi)];
end
L = legend(legend_matrix);
set(L,'Units','characters');
hold off

x_lim = xlim;
y_lim = ylim;
T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'scale_cv_functions',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% SPIKE TRIGGERED WAVELET
% --------------------------------------------------------------------
function spike_triggered_wavelet(h,handles,fdi_index,strinp,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end
    
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method '\stw\' strinp];
fln = ['STW_' method '_' CELL '_substate' num2str(fdi_index) '.jpg'];
ff = [pth '\' fln];
str = ['SpikeTriggeredWaveletImage = imread(''' ff ''');'];
eval(str)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Find the different states
fdi = find(diff(IntraBurstIvCv));
fdi = fdi + 1;
lenfdi = length(fdi);

% Create buttons
sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
if isequal(gcf,handles.figure1)
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''spike_triggered_wavelet'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw sca spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
else
    for i = 1:lenfdi
        spec = ['b_ica_gui4(''spike_triggered_wavelet'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
        button_struct(i).string = num2str(i);
        button_struct(i).callback = [sfw spec];
    end
    buts = create_buttons(h,handles,gca,button_struct);
end
if ~isempty(buts)
    fw = findobj(buts,'FontWeight','bold');
    if isempty(fw)
        set_font_weight(h,handles,gca,buts(1))
    end
end

% Plot
hold off
imagesc(SpikeTriggeredWaveletImage)
y_lim = ylim;
x_lim = xlim;

T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'spike_triggered_wavelet',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% SPIKE TRIGGERED WAVELET AVERAGE
% --------------------------------------------------------------------
function spike_triggered_wavelet_average(h,handles,strinp,varargin)

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
if ~isempty(varargin)       % get clustering method if called from open GUI
    method = varargin{1};
else                        % get clustering method if called from ICA GUI
    method = handles.clustering_method;
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_STW_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));
if strcmp(strinp,'power')
    StwMean = StwPowerMean;
    StwStd = StwPowerStd;
elseif strcmp(strinp,'phase')
    StwMean = StwPhaseMean;
    StwStd = StwPhaseStd;
end

global DATAPATH
ff = [DATAPATH 'ICA\ica_newave3\wavelet_settings\f.mat'];
load(ff)

% Store clustering method in axes application data
setappdata(gca,'clustering_method',method)

% Plot
lenfdi = length(StwMean);
hold off
clr = zeros(1,lenfdi);
for b = 1:lenfdi
    clr(b) = (b - 1) / (lenfdi - 1);
    semilogx(ScaleVector,StwMean{b},'Color',[0 1-clr(b) clr(b)]);
    hold on
%     semilogx(ScaleVector,StwMean{b}+StwStd{b},'Color','black');
%     semilogx(ScaleVector,StwMean{b}-StwStd{b},'Color','black');
    x_lim = xlim;
    y_lim = ylim;
    axis([1 ScaleVector(1) y_lim(1) y_lim(2)]);
    legend_matrix{b} = ['substate: ' num2str(b) '/' num2str(lenfdi)];
end
L = legend(legend_matrix);
set(L,'Units','characters');
hold off

x_lim = xlim;
y_lim = ylim;
T = text(x_lim(1)+4*(x_lim(2)-x_lim(1))/5,y_lim(1)+9*(y_lim(2)-y_lim(1))/10,...
    [method],'Color','k');      % Method
set(T,'Units','normalized');

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'spike_triggered_wavelet_average',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% SINGLE POINT SPIKE TRIGGERED WAVELET
% --------------------------------------------------------------------
function spstw(h,handles,long_index,strinp,varargin)

% Define long_index if called from OPEN GUI
if exist('long_index') == 0
    long_index = 1;
end

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\spstw\' strinp];
dd = dir(pth);
kk = zeros(1,length(dd));
for i = 3:length(dd)
    kk(i) = strcmp(dd(i).name(7:12),CELL(1:6));
end
nn = find(kk);
fln = dd(nn(long_index)).name;
needbuts = 0;
if length(nn) == 2
    needbuts = 1;
elseif length(nn) > 2
    error('Long registration.')
end
ff = [pth '\' fln];
str = ['SinglePointSpikeTriggeredWaveletImage = imread(''' ff ''');'];
eval(str)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Create buttons
if needbuts
    sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    if isequal(gcf,handles.figure1)
        for i = 1:2
            spec = ['b_ica_gui4(''spstw'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw sca spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    else
        for i = 1:2
            spec = ['b_ica_gui4(''spstw'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    end
    if ~isempty(buts)
        fw = findobj(buts,'FontWeight','bold');
        if isempty(fw)
            set_font_weight(h,handles,gca,buts(1))
        end
    end
end

% Plot
hold off
imagesc(SinglePointSpikeTriggeredWaveletImage)

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'spstw',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% SINGLE POINT SPIKE TRIGGERED WAVELET AVERAGE
% --------------------------------------------------------------------
function spstw_average(h,handles,long_index,strinp,varargin)

% Define long_index if called from OPEN GUI
if exist('long_index') == 0
    long_index = 1;
end

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\spstw'];
dd = dir(pth);
kk = zeros(1,length(dd));
for i = 3:length(dd)
    if ~dd(i).isdir
        kk(i) = strcmp(dd(i).name(17:22),CELL(1:6));
    else
        kk(i) = logical(0);
    end
end
nn = find(kk);
fln = dd(nn(long_index)).name;
needbuts = 0;
if length(nn) == 2
    needbuts = 1;
elseif length(nn) > 2
    error('Long registration.')
end
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

if strcmp(strinp,'power')
    SpstwMean = SpstwPowerMean;
    SpstwStd = SpstwPowerStd;
elseif strcmp(strinp,'phase')
    SpstwMean = SpstwPhaseMean;
    SpstwStd = SpstwPhaseStd;
end

global DATAPATH
ff = [DATAPATH 'ICA\ica_newave3\wavelet_settings\f.mat'];
load(ff)

% Create buttons
if needbuts
    sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    if isequal(gcf,handles.figure1)
        for i = 1:2
            spec = ['b_ica_gui4(''spstw_average'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw sca spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    else
        for i = 1:2
            spec = ['b_ica_gui4(''spstw_average'',gcbo,guidata(gcbo),' num2str(i) ',''' strinp ''')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    end
    if ~isempty(buts)
        fw = findobj(buts,'FontWeight','bold');
        if isempty(fw)
            set_font_weight(h,handles,gca,buts(1))
        end
    end
end

% Plot
hold off
semilogx(ScaleVector,SpstwMean,'Color','red');
% hold on
% semilogx(ScaleVector,SpstwMean+SpstwStd,'Color','black');
% semilogx(ScaleVector,SpstwMean-SpstwStd,'Color','black');
% hold off
x_lim = xlim;
y_lim = ylim;
axis([1 ScaleVector(1) y_lim(1) y_lim(2)]);

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'spstw_average',strinp)
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% THRESHOLDING
% --------------------------------------------------------------------
function thresholding(h,handles,long_index,varargin)

% Define long_index if called from OPEN GUI
if exist('long_index') == 0
    long_index = 1;
end

% Get 'handles' structure if called from undocked figure
if isempty(handles)
    global ICA_GUI_HANDLE
    h = ICA_GUI_HANDLE;
    handles = guidata(h);
    isundocked = 1;
else
    isundocked = 0;
end

% Load data
global CELL
if isempty(CELL)
    warndlg('No cell has been choosen. Ctrl + O for Open.','Warning');
    return
end

global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\thres'];
dd = dir(pth);
kk = zeros(1,length(dd));
for i = 3:length(dd)
    if ~dd(i).isdir
        kk(i) = strcmp(dd(i).name(17:22),CELL(1:6));
    else
        kk(i) = logical(0);
    end
end
nn = find(kk);
fln = dd(nn(long_index)).name;
needbuts = 0;
if length(nn) == 2
    needbuts = 1;
elseif length(nn) > 2
    error('Long registration.')
end
ff = [pth '\' fln];
load(ff)
fs = findstr(fln,'_');
fname = fln(fs(3)+1:fs(5)-1);
datinx1 = str2num(fln(fs(5)+1:fs(6)-1));
datinx2 = str2num(fln(fs(6)+1:end-4));

global DATADIR      % we have to load raw data: it slows down the display!
pth = [DATADIR];
dd = dir(pth);
kk = zeros(1,length(dd));
for i = 3:length(dd)
    if ~dd(i).isdir
        kk(i) = strcmp(dd(i).name(1:6),CELL(1:6));
    else
        kk(i) = logical(0);
    end
end
nn = find(kk);
fln = dd(nn).name;
ff = [pth '\' fln];
data = load(ff);
if isstruct(data)
    field = fieldnames(data);
    s = struct('type','.','subs',field);
    data = subsref(data,s);
end
if size(data,2) == 1
    data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)]; 
end
unit = data(datinx1:datinx2,2)';

% Create buttons
if needbuts
    sfw = ['b_ica_gui4(''set_font_weight'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    sca = ['b_ica_gui4(''set_current_axes'',gcbo,guidata(gcbo),' num2str(gca) ');'];
    if isequal(gcf,handles.figure1)
        for i = 1:2
            spec = ['b_ica_gui4(''thresholding'',gcbo,guidata(gcbo),' num2str(i) ')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw sca spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    else
        for i = 1:2
            spec = ['b_ica_gui4(''thresholding'',gcbo,guidata(gcbo),' num2str(i) ')'];
            button_struct(i).string = num2str(i);
            button_struct(i).callback = [sfw spec];
        end
        buts = create_buttons(h,handles,gca,button_struct);
    end
    if ~isempty(buts)
        fw = findobj(buts,'FontWeight','bold');
        if isempty(fw)
            set_font_weight(h,handles,gca,buts(1))
        end
    end
end

% Create 'time' vector
dt = 0.0001;
len = datinx2 - datinx1 + 1;
time = [1:len] * dt;

% Create pseudounit
if datinx1 > 2000000
    pseudoindex = datinx1 - 2000000;
    vdisc = vdisc - pseudoindex;
else
    vdisc = vdisc - datinx1;
end

z = zeros(1,length(time));
z(vdisc) = 1;
z(vdisc+1) = -1;

% Plot
hold off
plot(time,unit,'b');
y_lim = ylim;
axis([time(1) time(end) y_lim(1) y_lim(2)])
hold on

ind2 = 0;
next = 1;
while ind2 < len
    ind1 = ind2 + 1;
    ind2 = ind2 + seglen;
    
    line([time(ind1) time(min(ind2,len))],[T(next) T(next)],'Color','g')    % plot threshold
    next = next + 1;
end

% Set axes properties
set_axes_properties(h,handles)

% Set 'display' field
set_display_field(h,handles,'thresholding')
handles = guidata(h);

% Set keypress function if undocked
if isundocked
    set_keypressfcn(h,handles)
    set_zoom(h,handles)
end



% --------------------------------------------------------------------
% CLUSTERING menu callbacks
% --------------------------------------------------------------------
function varargout = single_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
set(handles.single_submenu,'Checked','on')
set(handles.average_submenu,'Checked','off')
set(handles.ward_submenu,'Checked','off')

% Set 'custering_method' field
handles.clustering_method = 'single';
guidata(h,handles)

% --------------------------------------------------------------------
function varargout = average_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
set(handles.single_submenu,'Checked','off')
set(handles.average_submenu,'Checked','on')
set(handles.ward_submenu,'Checked','off')

% Set 'custering_method' field
handles.clustering_method = 'average';
guidata(h,handles)

% --------------------------------------------------------------------
function varargout = ward_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
set(handles.single_submenu,'Checked','off')
set(handles.average_submenu,'Checked','off')
set(handles.ward_submenu,'Checked','on')

% Set 'custering_method' field
handles.clustering_method = 'ward';
guidata(h,handles)

% --------------------------------------------------------------------
function varargout = inter_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.inter_submenu)

% --------------------------------------------------------------------
function varargout = intra_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.intra_submenu)

% --------------------------------------------------------------------
function varargout = extra_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.extra_submenu)

% --------------------------------------------------------------------
function varargout = first_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.first_submenu)

% --------------------------------------------------------------------
function varargout = allfirst_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.allfirst_submenu)

% --------------------------------------------------------------------
function varargout = hold_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.hold_submenu)

% --------------------------------------------------------------------
function varargout = legend1_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.legend1_submenu)



% --------------------------------------------------------------------
% WAVELET AVERAGE menu callbacks
% --------------------------------------------------------------------
function varargout = legend2_submenu_Callback(h, eventdata, handles, varargin)

% Set 'Checked' property
checkswitch(handles.legend2_submenu)

% --------------------------------------------------------------------
function varargout = hz1_3_submenu_Callback(h, eventdata, handles, varargin)

% Set 'Checked' property
checkswitch(handles.hz1_3_submenu)

% --------------------------------------------------------------------
function varargout = hz3_6_submenu_Callback(h, eventdata, handles, varargin)

% Set 'Checked' property
checkswitch(handles.hz3_6_submenu)

% --------------------------------------------------------------------
function varargout = hz6_20_submenu_Callback(h, eventdata, handles, varargin)

% Set 'Checked' property
checkswitch(handles.hz6_20_submenu)

% --------------------------------------------------------------------
function varargout = hz20_50_submenu_Callback(h, eventdata, handles, varargin)

% Set 'Checked' property
checkswitch(handles.hz20_50_submenu)



% --------------------------------------------------------------------
% STRUCTURAL ELEMENTS menu callbacks
% --------------------------------------------------------------------
function varargout = superimpose1_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.superimpose1_submenu)

% --------------------------------------------------------------------
function varargout = superimpose2_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.superimpose2_submenu)

% --------------------------------------------------------------------
function varargout = superimpose3_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.superimpose3_submenu)

% --------------------------------------------------------------------
function varargout = superimpose4_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.superimpose4_submenu)

% --------------------------------------------------------------------
function varargout = crossline_submenu_Callback(h,eventdata,handles,varargin)

% Set 'Checked' property
checkswitch(handles.crossline_submenu)



% --------------------------------------------------------------------
% CHECKSWITCH
% --------------------------------------------------------------------
function checkswitch(menuhandle)

% Set 'Checked' property
ischk = get(menuhandle,'Checked');
if strcmp(ischk,'on')
    set(menuhandle,'Checked','off')
else
    set(menuhandle,'Checked','on')
end



% --------------------------------------------------------------------
% KEYPRESS FUNCTION for UNDOCKED STRUCTURAL ELEMENTS - superimpose the 
% wavelet average vectors, crossline the bursts
% --------------------------------------------------------------------
function callgreyplot

% Calling the appropriate functions
inp = get(gcf,'CurrentCharacter');
switch inp
case 'y'    % y - plots wavelet average magnitudes from 1 to 3 hz against time on the spike train
    greyplot('r')
case 'x'    % x - plots wavelet average magnitudes from 3 to 6 hz against time on the spike train
    greyplot('b')
case 'c'    % c - plots wavelet average magnitudes from 6 to 20 hz against time on the spike train
    greyplot('g')
case 'v'    % v - plots wavelet average magnitudes from 20 to 50 hz against time on the spike train
    greyplot('m')
case 'l'    % l - draws one line per burst connecting the intraburst action potentials
    burstliner
end

% --------------------------------------------------------------------
function greyplot(s)

% Get handles structure
global ICA_GUI_HANDLE
h = ICA_GUI_HANDLE;
handles = guidata(h);

% Load wavelet average vectors
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\wavevectors\'];
fln = ['WAVEVECTOR_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

sw = size(Wavevec,2);
ttm = linspace(i_first/10000,i_second/10000,sw);

% Plot
hold on
switch s
case 'r'
    if isempty(handles.greyplot1_handle)
        wm = Wavevec(4,:) - mean(Wavevec(4,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot1_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot1_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot1_handle,'Visible','off')
        else
            set(handles.greyplot1_handle,'Visible','on')
            if ~isempty(handles.greyplot2_handle)
                set(handles.greyplot2_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot3_handle)
                set(handles.greyplot3_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot4_handle)
                set(handles.greyplot4_handle,'Visible','off')
            end
        end
    end
case 'b'
    if isempty(handles.greyplot2_handle)
        wm = Wavevec(3,:) - mean(Wavevec(3,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot2_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot2_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot2_handle,'Visible','off')
        else
            set(handles.greyplot2_handle,'Visible','on')
            if ~isempty(handles.greyplot1_handle)
                set(handles.greyplot1_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot3_handle)
                set(handles.greyplot3_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot4_handle)
                set(handles.greyplot4_handle,'Visible','off')
            end
        end
    end
case 'g'
    if isempty(handles.greyplot3_handle)
        wm = Wavevec(2,:) - mean(Wavevec(2,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot3_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot3_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot3_handle,'Visible','off')
        else
            set(handles.greyplot3_handle,'Visible','on')
            if ~isempty(handles.greyplot1_handle)
                set(handles.greyplot1_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot2_handle)
                set(handles.greyplot2_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot4_handle)
                set(handles.greyplot4_handle,'Visible','off')
            end
        end
    end
case 'm'
    if isempty(handles.greyplot4_handle)
        wm = Wavevec(1,:) - mean(Wavevec(1,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot4_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot4_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot4_handle,'Visible','off')
        else
            set(handles.greyplot4_handle,'Visible','on')
            if ~isempty(handles.greyplot1_handle)
                set(handles.greyplot1_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot2_handle)
                set(handles.greyplot2_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot3_handle)
                set(handles.greyplot3_handle,'Visible','off')
            end
        end
    end
end
hold off
guidata(h,handles)

% --------------------------------------------------------------------
function burstliner

% Get handles structure
global ICA_GUI_HANDLE
h = ICA_GUI_HANDLE;
handles = guidata(h);
method = getappdata(gca,'clustering_method');

% Load data
global CELL
global DATAPATH
pth = [DATAPATH 'ICA\ica_newave3\' method];
fln = ['THETA_ICA_UNIT_' method '_' CELL '.mat'];
ff = [pth '\' fln];
load(ff)
fs2 = findstr(CELL,'_');
i_first = str2num(CELL(fs2(2)+1:fs2(3)-1));
i_second = str2num(CELL(fs2(3)+1:end));

% Create 'time' vector
dt = 0.0001;
len = i_second - i_first;
time = [1:i_second] * dt; 

% Get pushbutton handles
ch = get(gcf,'Children');
button_handles = findobj(ch,'Style','pushbutton');

% Get current substate ('fdi_index')
ax = gca;
ax = round(ax*10000)/10000;
next = 0;
for i = 1:length(button_handles)
    usd = get(button_handles(i),'UserData');
    axes_handle = usd.axes_handle;
    axes_handle = round(axes_handle*10000)/10000;
    if isequal(axes_handle,ax)
        fw = get(button_handles(i),'FontWeight');
        if strcmp(fw,'bold')
            fdi_index = next;
            break
        else
            next = next + 1;
        end
    end
end

% Find the different states
fdi = find(diff(IntraBurstIvCv));
fdi = fdi + 1;
lenfdi = length(fdi);
dec = fdi(fdi_index);

% Drawing the lines or switching their visibility
hold on;
sbd = size(Burst{dec},2);
if isempty(handles.line_handles)
    line_handles = zeros(1,sbd);
    for j = 1:sbd
        rajz = [time(Vdisc(Burst{dec}(1,j))) time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
        line_handles(j) = line(rajz(1,:),rajz(2,:),'Color','k');
    end
    handles.line_handles = line_handles;
else
    vis = get(handles.line_handles,'Visible');
    if strcmp(vis{1},'on')
        set(handles.line_handles,'Visible','off')
    else
        set(handles.line_handles,'Visible','on')
    end
end
hold off
guidata(h,handles)