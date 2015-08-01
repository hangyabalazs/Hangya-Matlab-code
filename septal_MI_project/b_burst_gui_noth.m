function varargout = b_burst_gui_noth(varargin)
%B_BURST_GUI_NOTH   GUI for cluster analysis.
%   BURST_GUI_NOTH is a tool for cluster analysis on multiple data files.
%   For every datafile it is possible to choose an optimal cluster number
%   (from 2 to 5) - 0 stands for non-bursty cells. Result is saved in
%   results' directory. Edit the program code to modify input and result
%   directories! BURST_GUI_NOTH works on the longest non-theta segment!
%
%   The program is implemented for 3-channel data files as well.
%
%   BURST_GUI_NOTH runs on non-theta segments (theta analysis is required)!
%
%   Ward's hierarchic clustering method is used. See ITCLUST and BURSTRUN 
%   for details.
%
%   Left axes: return plot. Right side axes: 'structural elements' (press
%   'c' to crossline bursts; 'z' to switch zoom on). User is able to choose
%   between 'All' or 'New' options using the popupmenu: only previouly
%   unsaved data are shown when choosing 'New'. Switch between cells using
%   the buttons with arrows! Cluster number selection panel: type optimal
%   cluster number in the edit text box, then press 'Save'/'Overwrite'.
%
%   See also ITCLUST and BURSTRUN.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @b_burst_gui_noth_OpeningFcn, ...
                   'gui_OutputFcn',  @b_burst_gui_noth_OutputFcn, ...
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
% Executes just before b_burst_gui_noth is made visible.
% -------------------------------------------------------------------------
function b_burst_gui_noth_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for b_burst_gui_noth
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set pushbutton 'CData'
matlabroot1 = 'C:\MATLAB_R2007a';
Cstart = imread([matlabroot1 '\work\Balazs\doubleleft_arrow.bmp']);
Cprev = imread([matlabroot1 '\work\Balazs\left_arrow.bmp']);
Cnext = imread([matlabroot1 '\work\Balazs\right_arrow.bmp']);
Clast = imread([matlabroot1 '\work\Balazs\doubleright_arrow.bmp']);
set(handles.start,'CData',Cstart);
set(handles.next,'CData',Cnext);
set(handles.previous,'CData',Cprev);
set(handles.last,'CData',Clast);

% Initialize
global DATAPATH
handles.inpdir = [DATAPATH 'Burst\Burst\noth_hc\'];   % input directory
handles.resdir = [DATAPATH 'Burst\Cluster\Noth_hc\'];   % result directory
handles.cell_no = 1;        % first cell in input directory
handles.inxprop = 'real';       % indexing problem
guidata(hObject,handles)

% Store default GUI position in the handles structure for the Resize Function
old_units = get(handles.figure1,'Units');
set(handles.figure1,'Units','pixels');
default_gui_pos = get(handles.figure1,'Position');
set(handles.figure1,'Units',old_units);
handles.gui_pos = default_gui_pos;
guidata(hObject,handles);

% Fill axes
fill_axes(hObject,eventdata,handles)



% -------------------------------------------------------------------------
% Outputs from this function are returned to the command line.
% -------------------------------------------------------------------------
function varargout = b_burst_gui_noth_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;



% --------------------------------------------------------------------
% RESIZE function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
gui_pos = handles.gui_pos;

% Remember the old Units properties
old_units_fig = get(handles.figure1,'Units');
old_units_axes = get(handles.axes1,'Units');
old_units_edit = get(handles.edit1,'Units');
old_units_text = get(handles.text1,'Units');
old_units_push = get(handles.save,'Units');
old_units_popupmenu = get(handles.popupmenu1,'Units');
old_units_uipanel = get(handles.uipanel1,'Units');

% Set Units to pixels
set(handles.figure1,'Units','pixels');
pos = get(handles.figure1,'Position');

% Compensate if too narrow
nlim = 900;
if pos(3) < nlim
    if pos(1) == gui_pos(1)
        pos(3) = nlim;
        set(handles.figure1,'Position',pos);
    else
        pos(1) = pos(1) + pos(3) - nlim;
        pos(3) = nlim;
        set(handles.figure1,'Position',pos);
    end
end

% Compensate if too low
llim = 650;
if pos(4) < llim
    if pos(2) == gui_pos(2)
        pos(4) = llim;
        set(handles.figure1,'Position',pos);
    else
        pos(2) = pos(2) + pos (4) - llim;
        pos(4) = llim;
        set(handles.figure1,'Position',pos);
    end
end

% New position of the 'save' pushbutton
% set(handles.save,'Units','pixels');
% pos_push = get(handles.save,'Position');
% set(handles.save,'Position',[pos_push(1) pos_push(2)+pos(4)-gui_pos(4)...
%         pos_push(3) pos_push(4)]);
% set(handles.save,'Units',old_units_push);

% New position of the 'shift cell' pushbuttons
set(handles.start,'Units','pixels');
pos_push = get(handles.start,'Position');
set(handles.start,'Position',[pos_push(1) pos_push(2)+pos(4)-gui_pos(4)...
        pos_push(3) pos_push(4)]);
set(handles.start,'Units',old_units_push);

set(handles.previous,'Units','pixels');
pos_push = get(handles.previous,'Position');
set(handles.previous,'Position',[pos_push(1) pos_push(2)+pos(4)-gui_pos(4)...
        pos_push(3) pos_push(4)]);
set(handles.previous,'Units',old_units_push);

set(handles.next,'Units','pixels');
pos_push = get(handles.next,'Position');
set(handles.next,'Position',[pos_push(1) pos_push(2)+pos(4)-gui_pos(4)...
        pos_push(3) pos_push(4)]);
set(handles.next,'Units',old_units_push);

set(handles.last,'Units','pixels');
pos_push = get(handles.last,'Position');
set(handles.last,'Position',[pos_push(1) pos_push(2)+pos(4)-gui_pos(4)...
        pos_push(3) pos_push(4)]);
set(handles.last,'Units',old_units_push);

% New position of the title edit text
set(handles.edit1,'Units','pixels');
pos_edit = get(handles.edit1,'Position');
set(handles.edit1,'Position',[pos_edit(1)+(pos(3)-gui_pos(3))/2 pos_edit(2)+pos(4)-gui_pos(4) ...
        pos_edit(3) pos_edit(4)]);
set(handles.edit1,'Units',old_units_push);

% New position of the edit text box in the 'Cluster number selection' panel
% set(handles.optno,'Units','pixels');
% pos_edit = get(handles.optno,'Position');
% set(handles.optno,'Position',[pos_edit(1) pos_edit(2)+pos(4)-gui_pos(4)...
%         pos_edit(3) pos_edit(4)]);
% set(handles.optno,'Units',old_units_push);

% New position of the static text in the 'Cluster number selection' panel
% set(handles.text1,'Units','pixels');
% pos_text = get(handles.text1,'Position');
% set(handles.text1,'Position',[pos_text(1) pos_text(2)+pos(4)-gui_pos(4)...
%     pos_text(3) pos_text(4)]);
% set(handles.text1,'Units',old_units_text);

% New position of the static text boxes on the right side
set(handles.text2,'Units','pixels');
pos_text = get(handles.text2,'Position');
set(handles.text2,'Position',[pos_text(1)+pos(3)-gui_pos(3) pos_text(2)+3*(pos(4)-gui_pos(4))/4 ...
    pos_text(3) pos_text(4)]);
set(handles.text2,'Units',old_units_text);

set(handles.text3,'Units','pixels');
pos_text = get(handles.text3,'Position');
set(handles.text3,'Position',[pos_text(1)+pos(3)-gui_pos(3) pos_text(2)+2*(pos(4)-gui_pos(4))/4 ...
    pos_text(3) pos_text(4)]);
set(handles.text3,'Units',old_units_text);

set(handles.text4,'Units','pixels');
pos_text = get(handles.text4,'Position');
set(handles.text4,'Position',[pos_text(1)+pos(3)-gui_pos(3) pos_text(2)+(pos(4)-gui_pos(4))/4 ...
    pos_text(3) pos_text(4)]);
set(handles.text4,'Units',old_units_text);

set(handles.text5,'Units','pixels');
pos_text = get(handles.text5,'Position');
set(handles.text5,'Position',[pos_text(1)+pos(3)-gui_pos(3) pos_text(2)...
    pos_text(3) pos_text(4)]);
set(handles.text5,'Units',old_units_text);

% New position of the popup menu
set(handles.popupmenu1,'Units','pixels');
pos_popupmenu = get(handles.popupmenu1,'Position');
set(handles.popupmenu1,'Position',[pos_popupmenu(1) pos_popupmenu(2)+pos(4)-gui_pos(4)...
    pos_popupmenu(3) pos_popupmenu(4)]);
set(handles.popupmenu1,'Units',old_units_popupmenu);

% New position of the uipanel
set(handles.uipanel1,'Units','pixels');
pos_uipanel = get(handles.uipanel1,'Position');
set(handles.uipanel1,'Position',[pos_uipanel(1) pos_uipanel(2)+pos(4)-gui_pos(4)...
    pos_uipanel(3) pos_uipanel(4)]);
set(handles.uipanel1,'Units',old_units_uipanel);

% New position of the 'return plot' axes
set(handles.axes1,'Units','pixels');
pos_axes = get(handles.axes1,'Position');
set(handles.axes1,'Position',[pos_axes(1) pos_axes(2)+pos(4)-gui_pos(4)...
    pos_axes(3) pos_axes(4)]);
set(handles.axes1,'Units',old_units_axes);

% New position of the 'structural elements' axes
set(handles.axes2,'Units','pixels');
pos_axes = get(handles.axes2,'Position');
set(handles.axes2,'Position',[pos_axes(1) pos_axes(2)+3*(pos(4)-gui_pos(4))/4 ...
    pos_axes(3)+pos(3)-gui_pos(3) pos_axes(4)+(pos(4)-gui_pos(4))/4]);
set(handles.axes2,'Units',old_units_axes);

set(handles.axes3,'Units','pixels');
pos_axes = get(handles.axes3,'Position');
set(handles.axes3,'Position',[pos_axes(1) pos_axes(2)+2*(pos(4)-gui_pos(4))/4 ...
    pos_axes(3)+pos(3)-gui_pos(3) pos_axes(4)+(pos(4)-gui_pos(4))/4]);
set(handles.axes3,'Units',old_units_axes);

set(handles.axes4,'Units','pixels');
pos_axes = get(handles.axes4,'Position');
set(handles.axes4,'Position',[pos_axes(1) pos_axes(2)+(pos(4)-gui_pos(4))/4 ...
    pos_axes(3)+pos(3)-gui_pos(3) pos_axes(4)+(pos(4)-gui_pos(4))/4]);
set(handles.axes4,'Units',old_units_axes);

set(handles.axes5,'Units','pixels');
pos_axes = get(handles.axes5,'Position');
set(handles.axes5,'Position',[pos_axes(1) pos_axes(2)...
    pos_axes(3)+pos(3)-gui_pos(3) pos_axes(4)+(pos(4)-gui_pos(4))/4]);
set(handles.axes5,'Units',old_units_axes);

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% -------------------------------------------------------------------------
% SAVE
% -------------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)

% Get optimal cluster number
ClusterNumber = str2num(get(handles.optno,'String'));

% Load
load(handles.loaded)

% Save
fn = [handles.resdir handles.noth_name '_CLUSTER'];
str1 = ['save ' fn ' IntraBurstIvCv ExtraBurstIvCv InterBurstIvCv'];
str2 = [' FirstSpikeCv AllFirstSpikeCv BurstLengthCv Burstiness'];
str3 = [' IntraBurstFrequency IntraBurstSpikeNumber BurstLength'];
str4 = [' BurstFrequency Burst Ccc ReturnPlotXData ReturnPlotYData'];
str5 = [' IdX Centroid SumDist'];
str6 = [' Vdisc'];
str7 = [' ClusterNumber'];
str = [str1 str2 str3 str4 str5 str6 str7];
eval(str)

% Set cluster no. selection background color
set(handles.optno,'BackGroundColor','c')

% Declare cell number as 'pseudo' - problem of indexing
if isequal(get(handles.popupmenu1,'Value'),2)    % if 'New' option is set
    handles.inxprop = 'pseudo';
    guidata(hObject,handles)
end

% -------------------------------------------------------------------------
function optno_Callback(hObject, eventdata, handles)



% -------------------------------------------------------------------------
% SHIFT CELL
% -------------------------------------------------------------------------
function previous_Callback(hObject,eventdata,handles)

% Shift
handles.cell_no = handles.cell_no - 1;
guidata(hObject,handles)
fill_axes(hObject,eventdata,handles)
handles = guidata(hObject);

% Declare cell number as 'real' - problem of indexing
handles.inxprop = 'real';
guidata(hObject,handles)

% -------------------------------------------------------------------------
function next_Callback(hObject,eventdata,handles)

% Shift
if isequal(handles.inxprop,'real')
    handles.cell_no = handles.cell_no + 1;
    guidata(hObject,handles)
end
fill_axes(hObject,eventdata,handles)
handles = guidata(hObject);

% Declare cell number as 'real' - problem of indexing
handles.inxprop = 'real';
guidata(hObject,handles)

% -------------------------------------------------------------------------
function start_Callback(hObject,eventdata,handles)

% Shift
handles.cell_no = 1;
guidata(hObject,handles)
fill_axes(hObject,eventdata,handles)
handles = guidata(hObject);

% Declare cell number as 'real' - problem of indexing
handles.inxprop = 'real';
guidata(hObject,handles)

% -------------------------------------------------------------------------
function last_Callback(hObject,eventdata,handles)

% Shift
handles.cell_no = 'last';
guidata(hObject,handles)
fill_axes(hObject,eventdata,handles)
handles = guidata(hObject);

% Declare cell number as 'real' - problem of indexing
handles.inxprop = 'real';
guidata(hObject,handles)



% -------------------------------------------------------------------------
% FILL AXES
% -------------------------------------------------------------------------
function fill_axes(hObject,eventdata,handles)

% Find longest theta and nontheta segment
files = dir(handles.inpdir);    % creat file lists
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
fn_long = {};
fn_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        fn_long{end+1} = files2(end).name;
        fn_short{end+1} = files2(end).name(1:6);
    end
end
files2 = files2(2:end);
sf = length(files2);

liststring = unique(fn_short);     % create liststring (cells in input dir.)
k2 = dir(handles.resdir);    % create readystring (cells in output dir.)
lenlist2 = length(k2)-2;
readystring = {};
readystring_long = {};
for t = 1:lenlist2
    if ~k2(t+2).isdir
        readystring{end+1} = k2(t+2).name(1:6);
        readystring_long{end+1} = k2(t+2).name;
    end
end
if isequal(get(handles.popupmenu1,'Value'),2)    % if 'New' option is set
    cell_list = setdiff(liststring,readystring);
    if isempty(cell_list)
        msgbox('All files in input directory has been processed.','Ready')
        return
    end
else
    cell_list = liststring;
end

if isequal(handles.cell_no,'last')      % cell name
    cname = cell_list{end};
else
    if handles.cell_no < 1
        handles.cell_no = 1;
    elseif handles.cell_no > length(cell_list)
        handles.cell_no = length(cell_list);
    end
    cname = cell_list{handles.cell_no};
end

if ismember(cname,readystring)      % all or new, overwrite
    set(handles.optno,'BackGroundColor','y')
    inx = find(strcmp(cname,readystring));
    rname = readystring_long{inx};
    str = ['load(''' handles.resdir rname ''')'];
    eval(str);
    set(handles.optno,'String',ClusterNumber);
    set(handles.save,'String','Overwrite')
else
    set(handles.optno,'BackGroundColor','w')
    set(handles.optno,'String','');
    set(handles.save,'String','Save')
end

fnd = find(strcmp(cname,fn_short));       % theta and non-theta list
act_list = fn_long(fnd);
noth_list = struct([]);
for i = 1:length(act_list)
    if ~isempty(findstr('noth',act_list{i}))
        noth_list(end+1).name = act_list{i};
        cmps = strread(noth_list(end).name,'%s','delimiter','_MH'); % 3ch data needs M and H delimiter
        noth_list(end).length = str2num(cmps{5}) - str2num(cmps{4});
        noth_list(end).i_first = str2num(cmps{4});
        noth_list(end).i_second = str2num(cmps{5});
    end
end
len_list = [noth_list.length];
mxlc = find(len_list==max(len_list));
noth = noth_list(mxlc).name;
cmps2 = strread(noth,'%s','delimiter','_');
handles.noth_name = cmps2{1};
for nn = 2:length(cmps2)-1
    handles.noth_name = [handles.noth_name '_' cmps2{nn}];
end

% Load
handles.loaded = [handles.inpdir noth];
load(handles.loaded)
guidata(hObject,handles)

% Create 'time' vector
dt = 0.0001;
time = [1:noth_list(mxlc).i_second] * dt;

% Plot
axes(handles.axes1)
plot(ReturnPlotXData,ReturnPlotYData,'b.','MarkerSize',8)

for dec = 2:5
    axes(eval(['handles.axes' num2str(dec)]))
    z = zeros(1,noth_list(mxlc).i_second);
    z(Vdisc) = 1;
    z(Vdisc+1) = -1;
    hold off
    plot_handle = plot(time(noth_list(mxlc).i_first:noth_list(mxlc).i_second),...
        z(noth_list(mxlc).i_first:noth_list(mxlc).i_second),'b');
    axis([time(noth_list(mxlc).i_first) time(noth_list(mxlc).i_second) -1.5 1.5])
    hold on;
    sbd = size(Burst{dec},2);
    hans = zeros(1,sbd);
    for j = 1:sbd
        ind1 = Vdisc(Burst{dec}(1,j));
        ind2 = Vdisc(Burst{dec}(2,j));
        hans(j) = plot(time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r');
    end
    hold off
    
    if isappdata(gca,'iscrossed')       % crossline if previous was crossed
        iscrossed = getappdata(gca,'iscrossed');
        if iscrossed
            line_handles = zeros(1,sbd);
            for j = 1:sbd
                rinx1 = time(Vdisc(Burst{dec}(1,j)));
                rinx2 = time(Vdisc(Burst{dec}(2,j)));
                rajz = [rinx1 rinx2;0.2 0.7];
                line_handles(j) = line(rajz(1,:),rajz(2,:),'Color','k');
            end
            setappdata(gca,'line_handles',line_handles);
        end
    end
    
    setappdata(gca,'burst',Burst{dec});      % set application data
    setappdata(gcf,'Vdisc',Vdisc);
    setappdata(gcf,'time',time);
    
    ach = get(gca,'Children');       % set buttondown function
    set([gca; ach],'ButtonDownFcn','b_burst_gui_noth(''ccax'',gcbo,[],guidata(gcbo))');
end

% Set keypress function
set(gcf,'KeypressFcn','b_burst_gui_noth(''KPFcn'',gcbo,[],guidata(gcbo))');
guidata(hObject,handles)



% -------------------------------------------------------------------------
% KEYPRESS function for axes
% -------------------------------------------------------------------------
function KPFcn(hObject,eventdata,handles)

inp = get(gcf,'CurrentCharacter');
switch inp
case 'c'
    crossline(hObject,eventdata,handles)
case 'z'
    zoom xon
end

% -------------------------------------------------------------------------
function crossline(h,eventdata,handles)

% Get application data
burst = getappdata(gca,'burst');
Vdisc = getappdata(gcf,'Vdisc');
time = getappdata(gcf,'time');

% Plot
hold on;        % crossline
sbd = size(burst,2);
if ~isappdata(gca,'line_handles')
    line_handles = zeros(1,sbd);
    for j = 1:sbd
        rajz = [time(Vdisc(burst(1,j))) time(Vdisc(burst(2,j)));0.2 0.7];
        line_handles(j) = line(rajz(1,:),rajz(2,:),'Color','k');
    end
    setappdata(gca,'line_handles',line_handles);
    setappdata(gca,'iscrossed',1)
else
    line_handles = getappdata(gca,'line_handles');
    vis = get(line_handles,'Visible');
    if strcmp(vis{1},'on')
        set(line_handles,'Visible','off')
        setappdata(gca,'iscrossed',0)
    else
        set(line_handles,'Visible','on')
        setappdata(gca,'iscrossed',1)
    end
end
hold off
guidata(h,handles)



% -------------------------------------------------------------------------
% BUTTONDOWN function for axes
% -------------------------------------------------------------------------
function ccax(hObject,eventdata,handles)

% Change current axes



% -------------------------------------------------------------------------
% POPUPMENU callback - 'All' or 'New'
% -------------------------------------------------------------------------
function popupmenu1_Callback(hObject, eventdata, handles)

% Jump to first cell
handles.cell_no = 1;
guidata(hObject,handles)

% Fill axes
fill_axes(hObject,eventdata,handles)



% -------------------------------------------------------------------------
% CREATE functions
% -------------------------------------------------------------------------
function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function optno_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end