function varargout = b_ica_gui2(varargin)
%ICA_GUI2   Graphical User Interface for Iterative Cluster Analysis.
%   ICA_GUI2 shows the result of variance analysis in the upper axes and return plot
%   or spike train in the lower axes.
%
%   Variance (ICA) results: variance values against number of clusters. Calculated 
%   variances are interburst interval, intraburst interval, extraburst interval,
%   inter-first-spike interval, inter-all-first-spike interval (first spikes and 
%   single spikes) variance. You can exclude each of them using the checkboxes, or
%   you can switch hold property of the upper axes with 'HOLD' checkbox.
%
%   Return plot or spike train: use the lower listbox or the mouse buttons to 
%   switch between return plot or spike train for different cluster numbers.
%   Intraburst spikes are red.
%
%   Use the upper listbox or the mouse buttons to switch between registration
%   segments.
%
%   You can set the input directory in 'Options' - 'Set input' menu. See ICA_ANAL
%   for details.
%
%   You are able to export figures from the GUI by pressing 'e' and 'r' buttons.
%   See CALLEXPORTFIG0 for details.
%
%   See also ICA_BETA2 and ICA_ANAL.

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
    
    % Set default value of setinput
    handles.setinput = 'base';
    
    % Set default value of varcell
    handles.varcell = [];
    
    eventdata = [];
    setliststring(fig, eventdata, handles, varargin);
    
    % Store default GUI position in the handles structure for the Resize Function
    old_units = get(handles.figure1,'Units');
    set(handles.figure1,'Units','pixels');
    default_gui_pos = get(handles.figure1,'Position');
    set(handles.figure1,'Units',old_units);
    handles.gui_pos = default_gui_pos;    
    guidata(fig,handles);
    
    % Set default value of checkboxes to 1
    set(handles.checkbox1,'Value',1)
    set(handles.checkbox2,'Value',1)
    set(handles.checkbox3,'Value',1)
    set(handles.checkbox4,'Value',1)
    set(handles.checkbox5,'Value',1)


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
% ALLVAR LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)
if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
    lb(h, eventdata, handles, varargin)
end

% --------------------------------------------------------------------
function varargout = lb(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');

fln = file_list{index_selected};
i_first = str2num(fln(11:16));
i_second = str2num(fln(18:23));

if isequal(handles.setinput,['base'])
    nm = '';
else
    nm = ['ica_beta_' handles.setinput];
end
ff = fullfile(DATAPATH,'ICA\ica_gui2',nm,file_list{index_selected});
load(ff)

% Plot on axis1 - variance plot
axes(handles.axes1)
bdf1 = get(handles.axes1,'ButtonDownFcn');
d = length(IntraBurstIvVar);

isintra = get(handles.checkbox1,'Value');
isextra = get(handles.checkbox2,'Value');
isinter = get(handles.checkbox3,'Value');
isfirst = get(handles.checkbox4,'Value');
isallfirst = get(handles.checkbox5,'Value');
isholdchecked = get(handles.checkbox6,'Value');

if isintra
    plot1_handle = plot([1:d],IntraBurstIvVar,'g');
    hold on
end
if isextra
    plot1_handle = plot([1:d],ExtraBurstIvVar,'b');
    hold on
end
if isinter
    plot1_handle = plot([1:d],InterBurstIvVar,'c');
    hold on
end
if isfirst
    plot1_handle = plot([1:d],FirstSpikeVar,'m');
    hold on
end
if isallfirst
    plot1_handle = plot([1:d],AllFirstSpikeVar,'k');
    hold on
end

legend('intraburstivvar','extraburstivvar','interburstivvar','firstspikevar','allfirstspikevar',2);
xlim([1 d]);
x_lim = xlim;
y_lim = ylim;
kk = (y_lim(2) - y_lim(1)) / 3;
length_vdisc = length(Vdisc);
if ~isholdchecked
    if length_vdisc < 100,
        text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length_vdisc)],'Color','r');
    else
        text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length_vdisc)],'Color','k');
    end    
    hold off
end

% Create list for listbox2
list = cell(1,d-1);
list{1} = 'Return Plot';
for t = 2:d-1
    list{t} = ['Number of clusters +1 = ',int2str(t+1)];
end
set(handles.listbox2,'String',list);

% Set edit string
switch handles.setinput
case {'under','over'}
    mfsp = min(FirstSpikeVar(3:8));
    set(handles.edit1,'String',mfsp)
case {'decrease','nodecrease'}
    set(handles.edit1,'String',handles.setinput)
case {'base'}
    set(handles.edit1,'String','')
case {'intraunder','intraover'}
    ibiv = IntraBurstIvVar(3);
    set(handles.edit1,'String',ibiv)
end        

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
index_selected2 = get(handles.listbox2,'Value');
item_list = get(handles.listbox1,'String');
if isequal(index_selected2,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    z(Vdisc) = 1;
    z(Vdisc+1) = -1;
    dec = index_selected2 + 1;
%     plot2_handle = plot(Time,z,'Color',[ 0.631 0.941 1.000 ]);
    plot2_handle = plot(Time,z,'r');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc(Burst{dec}(1,j));
        ind2 = Vdisc(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'b')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','k');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
handles.plot2_handle = plot2_handle;
set(handles.plot2_handle,'ButtonDownFcn',bdf2)
set(handles.axes2,'ButtonDownFcn',bdf2);

handles.plot1_handle = plot1_handle;
set(handles.plot1_handle,'ButtonDownFcn',bdf1)
set(handles.axes1,'ButtonDownFcn',bdf1);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig0')

% Store variables in handles structure
MaxClusterNumber = d;
VarCell = cell(1,11);
VarCell{1} = IntraBurstIvVar;
VarCell{2} = ExtraBurstIvVar;
VarCell{3} = InterBurstIvVar;
VarCell{4} = FirstSpikeVar;
VarCell{5} = AllFirstSpikeVar;
VarCell{6} = Vdisc;
VarCell{7} = ReturnPlotXData;
VarCell{8} = ReturnPlotYData;
VarCell{9} = Time;
VarCell{10} = Burst;
VarCell{11} = MaxClusterNumber;

handles.varcell = VarCell;
guidata(h,handles);


% --------------------------------------------------------------------
% SPIKE TRAIN - RETURN PLOT LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox2_Callback(h, eventdata, handles, varargin)
if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
    lb2(h, eventdata, handles, varargin)
end

% --------------------------------------------------------------------
function varargout = lb2(h, eventdata, handles, varargin)

% Load ICA data
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');

fln = file_list{index_selected};
i_first = str2num(fln(11:16));
i_second = str2num(fln(18:23));

if isempty(handles.varcell)
    if isequal(handles.setinput,['base'])
        nm = '';
    else
        nm = ['ica_beta_' handles.setinput];
    end
    global DATAPATH
    ff = fullfile(DATAPATH,'ICA\ica_gui2',nm,file_list{index_selected});
    load(ff)
else
    VarCell = handles.varcell;
    IntraBurstIvVar = VarCell{1};
    ExtraBurstIvVar = VarCell{2};
    InterBurstIvVar = VarCell{3};
    FirstSpikeVar = VarCell{4};
    AllFirstSpikeVar = VarCell{5};
    Vdisc = VarCell{6};
    ReturnPlotXData = VarCell{7};
    ReturnPlotYData = VarCell{8};
    Time = VarCell{9};
    Burst = VarCell{10};
    MaxClusterNumber = VarCell{11};
end

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
index_selected2 = get(handles.listbox2,'Value');
item_list = get(handles.listbox1,'String');
bdf2 = get(handles.axes2,'ButtonDownFcn');
if isequal(index_selected2,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    z(Vdisc) = 1;
    z(Vdisc+1) = -1;
    dec = index_selected2 + 1;
%     plot2_handle = plot(Time,z,'Color',[ 0.631 0.941 1.000 ]);
    plot2_handle = plot(Time,z,'r');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc(Burst{dec}(1,j));
        ind2 = Vdisc(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'b')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','k');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
handles.plot2_handle = plot2_handle;
set(handles.plot2_handle,'ButtonDownFcn',bdf2)
set(handles.axes2,'ButtonDownFcn',bdf2);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig0')

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
old_units_edit = get(handles.edit1,'Units');
old_units_check = get(handles.checkbox1,'Units');

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
pos_list = cell(1,2);
for t = 1:2
    eval(['set(handles.listbox',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_list{t} = get(handles.listbox',int2str(t),',''Position'');']);
    eval(['set(handles.listbox',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_list{t}(1) pos(4)-gui_pos(4)+pos_list{t}(2) pos_list{t}(3) pos_list{t}(4)]);']);
end

% New position of the checkbox
pos_check = cell(1,6);
for t = 1:6
    eval(['set(handles.checkbox',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_check{t} = get(handles.checkbox',int2str(t),',''Position'');']);
    eval(['set(handles.checkbox',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_check{t}(1) pos(4)-gui_pos(4)+pos_check{t}(2) pos_check{t}(3) pos_check{t}(4)]);']);
end

% New position of the edit text
set(handles.edit1,'Units','pixels');
pos_edit = get(handles.edit1,'Position');
set(handles.edit1,...
    'Position',[pos_edit(1) pos_edit(2) pos(3)-gui_pos(3)+pos_edit(3) pos_edit(4)]);

% New position of the static text above the listbox
pos_text = cell(1,3);
set(handles.text1,'Units','pixels');
pos_text{1} = get(handles.text1,'Position');
set(handles.text1,...
    'Position',[pos(3)-gui_pos(3)+pos_text{1}(1) pos(4)-gui_pos(4)+pos_text{1}(2) pos_text{1}(3) pos_text{1}(4)]);

% New position of the axes
pos_axes = cell(1,2);
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
set(handles.edit1,'Units',old_units_edit);
set(handles.checkbox1,'Units',old_units_check);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox5_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox6_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = options_menu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = setinput_submenu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
% BASE SUBMENU callback
% --------------------------------------------------------------------
function varargout = base_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'base';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



% --------------------------------------------------------------------
% UNDER SUBMENU callback
% --------------------------------------------------------------------
function varargout = under_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'under';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



% --------------------------------------------------------------------
% OVER SUBMENU callback
% --------------------------------------------------------------------
function varargout = over_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'over';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



% --------------------------------------------------------------------
% DECREASE SUBMENU callback
% --------------------------------------------------------------------
function varargout = decrease_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'decrease';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



% --------------------------------------------------------------------
% NODECREASE SUBMENU callback
% --------------------------------------------------------------------
function varargout = nodecrease_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'nodecrease';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



% --------------------------------------------------------------------
% INTRAUNDER SUBMENU callback
% --------------------------------------------------------------------
function varargout = intraunder_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'intraunder';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



% --------------------------------------------------------------------
% INTRAOVER SUBMENU callback
% --------------------------------------------------------------------
function varargout = intraover_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'intraover';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
guidata(h,handles)



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
if isequal(handles.setinput,['base'])
    nm = '';
else
    nm = ['ica_beta_' handles.setinput];
end
ff = fullfile(DATAPATH,'ICA\ica_gui2',nm);
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