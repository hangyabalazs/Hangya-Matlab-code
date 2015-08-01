function varargout = b_thres_gui(varargin)
%THRES_GUI    Graphical User Interface for thresholding and discrimination.
%   THRES GUI loads data preprocessed by THRESRUN2 (online version has not
%   been implemented yet). THRESRUN2 is an automatic thresholder using
%   THRES5. The unit with the multiple threshold is displayed in the axes.
%   User can zoom in an out with the mouse in the usual way. Pressing down
%   the shift key plus mouse click sets the threshold to new value. Same
%   result can be gained by editing the 'Threshold' box. Displayed segment
%   boundaries can be set int the 'Interval set' panel. Input directory
%   (THRESRUN2 output directory, saved threshold values) and output
%   directory (discriminated data files) can be set as well.
%   Accept button: accept displayed threshold, discriminate (with own disc
%   function) and save discriminated data.
%   Reject button: clear displayed threshold.
%   Thresholding failure rate (overall or for unique sessions) are displayed
%   on the GUI. Reset button sets 'This session' failure rate to zero.
%   The listbox displays all files in the input directory or files without
%   discriminated counterpair. User can switch between the above two modes
%   via the popup menu below the listbox. Online mode (without
%   preprocessing) has not been implemented yet.
%
%   See also THRESRUN2.

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
    
    % Directories
    global DATAPATH
    global DATADIR
    thresroot = fullfile(DATAPATH,'Thres\');
    thresconfig = fullfile(thresroot,'config\');
    try
        param = fullfile(thresconfig,'settings.mat');
        load(param)
        handles.inpdir = inpdir;
        handles.resdir = resdir;
    catch
        handles.inpdir = fullfile(DATAPATH,'Raphe\raphe_juxta\Thres\temp\');
        handles.resdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\temp\'];
        W1 = warndlg('Configuration file insufficient or not exist: input and result directory set to default.',...
            'Directory Configuration');
    end
    handles.datadir = [DATADIR 'PViktor\non-glycinergic\'];
    handles.sr = 10000;     % sampling rate
    set(handles.edit4,'String',handles.inpdir)
    set(handles.edit5,'String',handles.resdir)
    guidata(fig,handles);
    
    % Failure rate, All or new, Preprocessed or online
    try
        set(handles.popupmenu1,'Value',allornew)
        set(handles.popupmenu2,'Value',preoron)
        set(handles.popupmenu3,'Value',thisoroverall)
    catch
        W2 = warndlg('Popup menu configuration failed: configuration file insufficient or not exist.',...
            'Popup Menu Configuration');
    end
    
    % Initiate 'reject list' and 'all list'
    handles.reject_list = {};
    handles.all_list = {};
    try
        param = fullfile(thresconfig,'overall_failure_rate.mat');
        load(param)
        handles.overall_reject_list = overall_reject_list;
        handles.overall_all_list = overall_all_list;
    catch
        handles.overall_reject_list = {};
        handles.overall_all_list = {};
        W3 = warndlg('Configuration file insufficient or not exist: overall failure rate reset.',...
            'Failure Rate Configuration');
    end
    guidata(fig,handles)
    
    % Base of logarithm for zooming
    handles.la = 8;
    guidata(fig,handles)
    
    % Store default GUI position in the handles structure for the Resize Function
    old_units = get(handles.figure1,'Units');
    set(handles.figure1,'Units','pixels');
    default_gui_pos = get(handles.figure1,'Position');
    set(handles.figure1,'Units',old_units);
    handles.gui_pos = default_gui_pos;
    guidata(fig,handles);
    
    % List string
    eventdata = [];
    setliststring(fig,eventdata,handles,varargin);
    handles = guidata(fig);
    
    % Warning dialogs
    if exist('W1')
        figure(W1)
    end
    if exist('W2')
        figure(W2)
    end
    if exist('W3')
        figure(W3)
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
% RESIZE function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
gui_pos = handles.gui_pos;

% Remember the old Units properties
old_units_fig = get(handles.figure1,'Units');
old_units_axes = get(handles.axes1,'Units');
old_units_edit = get(handles.edit1,'Units');
old_units_text = get(handles.text1,'Units');
old_units_push = get(handles.pushbutton1,'Units');
old_units_listbox = get(handles.listbox1,'Units');
old_units_popupmenu = get(handles.popupmenu1,'Units');
old_units_frame = get(handles.frame1,'Units');

% Set Units to pixels
set(handles.figure1,'Units','pixels');
pos = get(handles.figure1,'Position');

% Compensate if too narrow
if pos(3) < 900
    if pos(1) == gui_pos(1)
        pos(3) = 900;
        set(handles.figure1,'Position',pos);
    else
        pos(1) = pos(1) + pos(3) - 900;
        pos(3) = 900;
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

% New position of the pushbuttons
pos_push = cell(1,3);
set(handles.pushbutton1,'Units','pixels');
pos_push{1} = get(handles.pushbutton1,'Position');
set(handles.pushbutton1,'Position',[pos_push{1}(1) pos_push{1}(2)...
        (pos(3)-gui_pos(3))/2+pos_push{1}(3) pos_push{1}(4)]);
set(handles.pushbutton1,'Units',old_units_push);

set(handles.pushbutton2,'Units','pixels');
pos_push{2} = get(handles.pushbutton2,'Position');
set(handles.pushbutton2,'Position',[(pos(3)-gui_pos(3))/2+pos_push{2}(1) pos_push{2}(2)...
        (pos(3)-gui_pos(3))/2+pos_push{2}(3) pos_push{2}(4)]);
set(handles.pushbutton2,'Units',old_units_push);

set(handles.pushbutton3,'Units','pixels');
pos_push{3} = get(handles.pushbutton3,'Position');
set(handles.pushbutton3,'Position',[(pos(3)-gui_pos(3))/3+pos_push{3}(1) pos_push{3}(2)...
        pos_push{3}(3) pos_push{3}(4)]);
set(handles.pushbutton3,'Units',old_units_push);

% New position of the listbox
set(handles.listbox1,'Units','pixels');
pos_list = get(handles.listbox1,'Position');
set(handles.listbox1,...
    'Position',[pos(3)-gui_pos(3)+pos_list(1) pos_list(2) pos_list(3) pos(4)-gui_pos(4)+pos_list(4)]);
set(handles.listbox1,'Units',old_units_listbox);

% New position of the edit text boxes in the 'Threshold' panel
pos_edit = cell(1,6);
for t = 1:3
    eval(['set(handles.edit',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_edit{t} = get(handles.edit',int2str(t),',''Position'');']);
    eval(['set(handles.edit',int2str(t),...
            ',''Position'',[(pos(3)-gui_pos(3))/3+pos_edit{t}(1) pos_edit{t}(2) pos_edit{t}(3) pos_edit{t}(4)]);']);
    eval(['set(handles.edit',int2str(t),',''Units'',old_units_edit);']);
end

% New position of the edit text boxes in the 'Directory' panel
for t = 4:5
    eval(['set(handles.edit',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_edit{t} = get(handles.edit',int2str(t),',''Position'');']);
    eval(['set(handles.edit',int2str(t),...
            ',''Position'',[2*(pos(3)-gui_pos(3))/3+pos_edit{t}(1) pos_edit{t}(2) pos_edit{t}(3) pos_edit{t}(4)]);']);
    eval(['set(handles.edit',int2str(t),',''Units'',old_units_edit);']);
end

% New position of the edit text in the 'Failure rate' panel
set(handles.edit6,'Units','pixels');
pos_edit{6} = get(handles.edit6,'Position');
set(handles.edit6,...
    'Position',[(pos(3)-gui_pos(3))/3+pos_edit{6}(1) pos_edit{6}(2) pos_edit{6}(3) pos_edit{6}(4)]);
set(handles.edit6,'Units',old_units_text);

% New position of the static text above the listbox
pos_text = cell(1,8);
set(handles.text1,'Units','pixels');
pos_text{1} = get(handles.text1,'Position');
set(handles.text1,...
    'Position',[pos(3)-gui_pos(3)+pos_text{1}(1) pos(4)-gui_pos(4)+pos_text{1}(2) pos_text{1}(3) pos_text{1}(4)]);
set(handles.text1,'Units',old_units_text);

% New position of the static text boxes in the 'Threshold' panel
for t = 2:5
    eval(['set(handles.text',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_text{t} = get(handles.text',int2str(t),',''Position'');']);
    eval(['set(handles.text',int2str(t),...
            ',''Position'',[(pos(3)-gui_pos(3))/3+pos_text{t}(1) pos_text{t}(2) pos_text{t}(3) pos_text{t}(4)]);']);
    eval(['set(handles.text',int2str(t),',''Units'',old_units_text);']);
end

% New position of the static text boxes in the 'Directory' panel
for t = 6:7
    eval(['set(handles.text',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_text{t} = get(handles.text',int2str(t),',''Position'');']);
    eval(['set(handles.text',int2str(t),...
            ',''Position'',[2*(pos(3)-gui_pos(3))/3+pos_text{t}(1) pos_text{t}(2) pos_text{t}(3) pos_text{t}(4)]);']);
    eval(['set(handles.text',int2str(t),',''Units'',old_units_text);']);
end

% New position of the static text in the 'Failure rate' panel
set(handles.text8,'Units','pixels');
pos_text{8} = get(handles.text8,'Position');
set(handles.text8,...
    'Position',[(pos(3)-gui_pos(3))/3+pos_text{8}(1) pos_text{8}(2) pos_text{8}(3) pos_text{8}(4)]);
set(handles.text8,'Units',old_units_text);

% New position of 'Threshold' panel frame
pos_frame = cell(1,3);
set(handles.frame1,'Units','pixels');
pos_frame{1} = get(handles.frame1,'Position');
set(handles.frame1,...
    'Position',[(pos(3)-gui_pos(3))/3+pos_frame{1}(1) pos_frame{1}(2) pos_frame{1}(3) pos_frame{1}(4)]);
set(handles.frame1,'Units',old_units_frame);

% New position of 'Directory' panel frame
set(handles.frame2,'Units','pixels');
pos_frame{2} = get(handles.frame2,'Position');
set(handles.frame2,...
    'Position',[2*(pos(3)-gui_pos(3))/3+pos_frame{2}(1) pos_frame{2}(2) pos_frame{2}(3) pos_frame{2}(4)]);
set(handles.frame2,'Units',old_units_frame);

% New position of 'Failure rate' panel frame
set(handles.frame3,'Units','pixels');
pos_frame{3} = get(handles.frame3,'Position');
set(handles.frame3,...
    'Position',[(pos(3)-gui_pos(3))/3+pos_frame{3}(1) pos_frame{3}(2) pos_frame{3}(3) pos_frame{3}(4)]);
set(handles.frame3,'Units',old_units_frame);

% New position of the popup menus under the listbox
pos_popupmenu = cell(1,3);
for t = 1:2
    eval(['set(handles.popupmenu',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_popupmenu{t} = get(handles.popupmenu',int2str(t),',''Position'');']);
    eval(['set(handles.popupmenu',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_popupmenu{t}(1) pos_popupmenu{t}(2) pos_popupmenu{t}(3) pos_popupmenu{t}(4)]);']);
    eval(['set(handles.popupmenu',int2str(t),',''Units'',old_units_popupmenu);']);
end

% New position of the popup menu in the 'Failure rate' panel
set(handles.popupmenu3,'Units','pixels');
pos_popupmenu{3} = get(handles.popupmenu3,'Position');
set(handles.popupmenu3,...
    'Position',[(pos(3)-gui_pos(3))/3+pos_popupmenu{3}(1) pos_popupmenu{3}(2) pos_popupmenu{3}(3) pos_popupmenu{3}(4)]);
set(handles.popupmenu3,'Units',old_units_popupmenu);

% New position of the axes
pos_axes = cell(1);
set(handles.axes1,'Units','pixels');
pos_axes{1} = get(handles.axes1,'Position');
set(handles.axes1,...
    'Position',[pos_axes{1}(1) pos_axes{1}(2) pos(3)-gui_pos(3)+pos_axes{1}(3)...
        (pos(4)-gui_pos(4))+pos_axes{1}(4)]);
set(handles.axes1,'Units',old_units_axes);

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% --------------------------------------------------------------------
% LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)

if isequal(get(handles.popupmenu2,'Value'),1)        % 'Preprocessed' mode
    lb_preprocessed(h,eventdata,handles,varargin)
elseif isequal(get(handles.popupmenu2,'Value'),2)    % 'On line' mode
    lb_online(h,eventdata,handles,varargin)
end

% --------------------------------------------------------------------
function varargout = lb_preprocessed(h, eventdata, handles, varargin)

% Load threshold
global DATAPATH
global DATADIR
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
fln = file_list{index_selected};
fn = ['THRES_' fln '.mat'];
ff = fullfile(handles.inpdir,fn);
load(ff)

% Load data
where2 = handles.datadir;
ff2 = [where2 fln];
try
    data0 = load(ff2);
    unitnames = {'Unit','unit'};
    eegnames = {'hEEG','EEG','eeg'};
    datafields = fieldnames(data0);
    eegname = intersect(eegnames,datafields);
    unitname = intersect(unitnames,datafields);
    unit = data0.(unitname{1}).values';
    eeg = data0.(eegname{1}).values';
catch
    data = b_load_data(ff2);
    unit = data(:,2)';
    eeg = data(:,1)';
end
handles.unit = unit;
handles.eeg = eeg;

min_max{1} = min(unit);
min_max{2} = max(unit);
min_max{3} = min(eeg);
min_max{4} = max(eeg);
handles.min_max = min_max;

handles.eegs = b_createvs(eeg,handles.la);
handles.units = b_createvs(unit,handles.la);

segfirst = 1;
seglast = length(unit);
handles.fname = fln;
handles.segfirst = segfirst;
handles.seglast = seglast;
handles.filename = [fln '.mat'];
guidata(h,handles)

% Set edit strings
set(handles.edit1,'String',num2str(segfirst/10000))
set(handles.edit2,'String',num2str(seglast/10000))
set(handles.edit3,'String','')

% Not rejected yet
handles.isrejected = 0;
guidata(h,handles)

% Create 'time' vector
dt = 0.0001;
len = seglast - segfirst + 1;
time = [1:len] * dt;
handles.time = time;
handles.len = len;
guidata(h,handles);

% Change GUI figure 'HandleVisibility'
default_handle_visibility = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','on');

% Plot unit
axes(handles.axes1)
handles.plot1_handle = plot(time,unit);
x_lim = xlim;
y_lim = ylim;
axis([time(1) time(end) y_lim(1) y_lim(2)]);
b_liner(handles.plot1_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);

% Plot threshold
ind2 = 0;
lenu = length(unit);
handles.threshold = [];
for i = 1:length(T)
    ind1 = ind2 + 1;
    ind2 = ind2 + seglen;
    handles.threshold(ind1:min(ind2,lenu)) = T(i);
end
guidata(h,handles)
plotthres(h,eventdata,handles,varargin)
handles = guidata(h);

% Set buttondown function
bdf1 = 'b_thres_gui(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
set(handles.plot1_handle,'ButtonDownFcn',bdf1)
set(handles.axes1,'ButtonDownFcn',bdf1);

% Reset GUI figure 'HandleVisibility'
set(handles.figure1,'HandleVisibility',default_handle_visibility);



% --------------------------------------------------------------------
function varargout = lb_online(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
% EDIT TEXT1 callback - set interval's LOWER LIMIT
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)

% Get input and limits
x1 = get(handles.edit1,'String');
x2 = get(handles.edit2,'String');
x_lim = get(handles.axes1,'XLim');
x1 = str2num(x1);
x2 = str2num(x2);
if isempty(x1)
    errordlg('Input must be numeric.','Error','modal');
    set(handles.edit1,'String',num2str(x_lim(1)));
    return
end
if x1 >= x2
    errordlg('Input values must be increasing.','Error','modal');
    set(handles.edit1,'String',num2str(x_lim(1)));
    return
end
y_lim = get(handles.axes1,'YLim');

% If out of range, set string back; else set axis and update plot
if x1 < handles.segfirst/10000
    set(handles.edit1,'String',num2str(handles.segfirst/10000));
    axes(handles.axes1)
    axis([handles.segfirst/10000 x_lim(2) y_lim(1) y_lim(2)]);
elseif x1 > handles.seglast/10000
    set(handles.edit1,'String',num2str(handles.seglast/10000));
    axes(handles.axes1)
    axis([handles.seglast/10000 x_lim(2) y_lim(1) y_lim(2)]);
else
    axes(handles.axes1)
    axis([x1 x_lim(2) y_lim(1) y_lim(2)]);
end
b_liner(handles.plot1_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
drawnow

% Clear threshold
set(handles.edit3,'String','');



% --------------------------------------------------------------------
% EDIT TEXT2 callback - set interval's UPPER LIMIT
% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)

% Get input and limits
x1 = get(handles.edit1,'String');
x2 = get(handles.edit2,'String');
x_lim = get(handles.axes1,'XLim');
x2 = str2num(x2);
x1 = str2num(x1);
if isempty(x2)
    errordlg('Input must be numeric.','Error','modal');
    set(handles.edit2,'String',num2str(x_lim(2)));
    return
end
if x1 >= x2
    errordlg('Input values must be increasing.','Error','modal');
    set(handles.edit2,'String',num2str(x_lim(2)));
    return
end
y_lim = get(handles.axes1,'YLim');

% If out of range, set string back; else set axis and update plot
if x2 < handles.segfirst/10000
    set(handles.edit2,'String',num2str(handles.segfirst/10000));
    axes(handles.axes1)
    axis([x_lim(1) handles.segfirst/10000 y_lim(1) y_lim(2)]);
elseif x2 > handles.seglast/10000
    set(handles.edit2,'String',num2str(handles.seglast/10000));
    axes(handles.axes1)
    axis([x_lim(1) handles.seglast/10000 y_lim(1) y_lim(2)]);
else
    axes(handles.axes1)
    axis([x_lim(1) x2 y_lim(1) y_lim(2)]);
end
b_liner(handles.plot1_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
drawnow

% Clear threshold
set(handles.edit3,'String','');



% --------------------------------------------------------------------
% EDIT TEXT3 callback - set THRESHOLD
% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)

% Get interval limits and threshold
x1 = get(handles.edit1,'String');
x2 = get(handles.edit2,'String');
if ~isempty(varargin)
    nth = varargin{1};   % call from axes button down function
else
    th = get(handles.edit3,'String');
    nth = str2num(th);
end

% Check input
if isempty(nth)
    errordlg('Input must be numeric.','Error','modal');
    set(handles.edit3,'String','');
    return
end

% Set threshold
nx1 = str2num(x1);
nx2 = str2num(x2);
% rnx1 = ceil(nx1);
% rnx2 = floor(nx2);
handles.threshold(nx1*10000:nx2*10000) = nth;
guidata(h,handles)

% Plot new threshold
plotthres(h,eventdata,handles,varargin)
handles = guidata(h);

% Rejected
handles.isrejected = 1;
guidata(h,handles)



% --------------------------------------------------------------------
% EDIT TEXT4 callback - INPUT DIRECTORY
% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
inpdir = get(handles.edit4,'String');
if ~isdir(inpdir)
    errordlg('Input directory do not exist.')
else
    handles.inpdir = inpdir;
end
guidata(h,handles);
setliststring(h,eventdata,handles,varargin);



% --------------------------------------------------------------------
% EDIT TEXT5 callback - RESULT DIRECTORY
% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)
resdir = get(handles.edit5,'String');
if ~isdir(resdir)
    errordlg('Result directory do not exist.')
else
    handles.resdir = resdir;
end
guidata(h,handles)
setliststring(h,eventdata,handles,varargin);



% --------------------------------------------------------------------
% EDIT TEXT6 callback - FAILURE RATE
% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
% PUSHBUTTON1 callback - ACCEPT
% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

% Check if sufficient threshold exists
if ~isempty(find(isnan(handles.threshold)))
    errordlg('Incomplete threshold','Threshold error')
    return
end

% Discrimination
vdisc = disc(handles.unit,handles.threshold);

% Get variables
eeg = handles.eeg;
fln = handles.filename;
resdir = handles.resdir;

% Save
fn = fullfile(resdir,[fln(1:end-4) '_d.mat']);
save(fn,'eeg','vdisc')

% Handle 'reject list' and failure rate
if handles.isrejected
    handles.reject_list{end+1} = fln;
    handles.reject_list = unique(handles.reject_list);
    handles.overall_reject_list{end+1} = fln;
    handles.overall_reject_list = unique(handles.overall_reject_list);
end
handles.all_list{end+1} = fln;
handles.all_list = unique(handles.all_list);
handles.overall_all_list{end+1} = fln;
handles.overall_all_list = unique(handles.overall_all_list);

pop = get(handles.popupmenu3,'Value');
if isequal(pop,1)
    rno = length(handles.reject_list);
    ano = length(handles.all_list);
    fr = rno / ano;
    set(handles.edit6,'String',num2str(fr))
elseif isequal(pop,2)
    rno = length(handles.overall_reject_list);
    ano = length(handles.overall_all_list);
    fr = rno / ano;
    set(handles.edit6,'String',num2str(fr))
end

% Plot next file
list = get(handles.listbox1,'String');
lenlist = length(list);
ix = get(handles.listbox1,'Value');
set(handles.listbox1,'Value',min(ix+1,lenlist))
listbox1_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
% PUSHBUTTON2 callback - REJECT
% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)

% Delete threshold from plot
if isfield(handles,'threshandle') & ishandle(handles.threshandle)
    delete(handles.threshandle)
end

% Delete threshold variable
if isfield(handles,'threshold')
    handles.threshold = repmat(NaN,1,length(handles.threshold));
    guidata(h,handles)
end

% Clear threshold edit text
set(handles.edit3,'String','');



% --------------------------------------------------------------------
% PUSHBUTTON3 callback - RESET FAILURE RATE
% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)

% Clear 'reject list' and 'all list'
handles.reject_list = {};
handles.all_list = {};
guidata(h,handles)

% Refresh failure rate
set(handles.edit6,'String','')



% --------------------------------------------------------------------
% POPUPMENU1 callback - SET LIST
% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)

% Set list
setliststring(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
% POPUPMENU3 callback - SET FAILURE RATE
% --------------------------------------------------------------------
function varargout = popupmenu3_Callback(h, eventdata, handles, varargin)

% Set failure rate
pop = get(handles.popupmenu3,'Value');
if isequal(pop,1)
    rno = length(handles.reject_list);
    ano = length(handles.all_list);
    fr = rno / ano;
    set(handles.edit6,'String',num2str(fr))
elseif isequal(pop,2)
    rno = length(handles.overall_reject_list);
    ano = length(handles.overall_all_list);
    fr = rno / ano;
    set(handles.edit6,'String',num2str(fr))
end



% --------------------------------------------------------------------
% Callback for Axes ButtonDownFunction - ZOOMING
% --------------------------------------------------------------------
function varargout = axes1_ButtonDownFcn(h, eventdata, handles, varargin)

% Get variables from 'handles' structure
time = handles.time;
min_max = handles.min_max;
min_unit = min_max{1};
max_unit = min_max{2};

% Set axis
seltyp = get(handles.figure1,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(handles.axes1,'CurrentPoint'); % button down detected
    units = get(handles.figure1,'units');
    set(handles.figure1,'units','pixels')
    rbbox([get(handles.figure1,'currentpoint') 0 0],get(handles.figure1,'currentpoint'),handles.figure1);                   % return figure units
    set(handles.figure1,'units',units)
    point2 = get(handles.axes1,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    if isequal(point1,point2)
        xx = get(handles.axes1,'XLim');
        yy = get(handles.axes1,'YLim');
        xx2 = (abs(xx(2) - xx(1))) / 4;
        yy2 = (abs(yy(2) - yy(1))) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        if time(1) < xx3(1) & time(end) > xx3(2)
            set(handles.axes1,'XLim',xx3);
        elseif time(1) > xx3(1)
            xx_new(1) = time(1);
            xx_new(2) = time(1) + (2 * xx2);
            set(handles.axes1,'XLim',xx_new);
        elseif time(end) < xx3(2)
            xx_new(1) = time(end) - (2 * xx2);
            xx_new(2) = time(end);
            set(handles.axes1,'XLim',xx_new);
        end
        if min_unit < yy3(1) & max_unit > yy3(2)
            set(handles.axes1,'YLim',yy3);
        elseif min_unit > yy3(1)
            yy_new(1) = min_unit;
            yy_new(2) = min_unit + (2 * yy2);
            set(handles.axes1,'YLim',yy_new);
        elseif max_unit < yy3(2)
            yy_new(1) = max_unit - (2 * yy2);
            yy_new(2) = max_unit;
            set(handles.axes1,'YLim',yy_new);
        end
    else
        axis([max(min([point1(1) point2(1)]),time(1)) min(max([point1(1) point2(1)]),time(end))...
                min([point1(2) point2(2)]) max([point1(2) point2(2)])]);
    end

    axes(handles.axes1);
    b_liner(handles.plot1_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
    
    
case 'open'   % set default
    axes(handles.axes1);
    axis([time(1) time(end) min_unit max_unit]);
    b_liner(handles.plot1_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
    
case 'extend'   % set threshold
    point = get(handles.axes1,'CurrentPoint'); % button down detected
    thres = point(1,2);
    edit3_Callback(h,eventdata,handles,thres)
    set(handles.edit3,'String',num2str(thres))
    
otherwise   % zoom out
    point = get(handles.axes1,'CurrentPoint'); % button down detected
    point = point(1,1:2);
    xx = get(handles.axes1,'XLim');
    yy = get(handles.axes1,'YLim');
    xx2 = abs(xx(2) - xx(1));
    yy2 = abs(yy(2) - yy(1));
    if xx2 > (time(end) - time(1)) / 2,
        xx2 = (time(end) - time(1)) / 2;
    end
    if yy2 > (abs(max_unit - min_unit)) / 2,
        yy2 = (abs(max_unit - min_unit)) / 2;
    end
    xx3(1) = point(1) - xx2;
    xx3(2) = point(1) + xx2;
    yy3(1) = point(2) - yy2;
    yy3(2) = point(2) + yy2;
    if time(1) < xx3(1) & time(end) > xx3(2)
        set(handles.axes1,'XLim',xx3);
    elseif time(1) > xx3(1)
        xx_new(1) = time(1);
        xx_new(2) = time(1) + (2 * xx2);
        set(handles.axes1,'XLim',xx_new);
    elseif time(end) < xx3(2)
        xx_new(1) = time(end) - (2 * xx2);
        xx_new(2) = time(end);
        set(handles.axes1,'XLim',xx_new);
    end
    if min_unit < yy3(1) & max_unit > yy3(2)
        set(handles.axes1,'YLim',yy3);
    elseif min_unit > yy3(1)
        yy_new(1) = min_unit;
        yy_new(2) = min_unit + (2 * yy2);
        set(handles.axes1,'YLim',yy_new);
    elseif max_unit < yy3(2)
        yy_new(1) = max_unit - (2 * yy2);
        yy_new(2) = max_unit;
        set(handles.axes1,'YLim',yy_new);
    end
    
    axes(handles.axes1);
    b_liner(handles.plot1_handle, handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
end
drawnow;

% Set edit strings
x_lim = get(handles.axes1,'XLim');
set(handles.edit1,'String',num2str(x_lim(1)));
set(handles.edit2,'String',num2str(x_lim(2)));



% --------------------------------------------------------------------
% Callback for GUI Close Request Function - SAVING PARAMETERS
% --------------------------------------------------------------------
function varargout = figure1_CloseRequestFcn(h, eventdata, handles, varargin)

try
    % Get variables
    inpdir = handles.inpdir;
    resdir = handles.resdir;
    allornew = get(handles.popupmenu1,'Value');
    preoron = get(handles.popupmenu2,'Value');
    thisoroverall = get(handles.popupmenu3,'Value');

    overall_reject_list = handles.overall_reject_list;
    overall_all_list = handles.overall_all_list;

    % Generate file names
    global DATAPATH
    thresroot = fullfile(DATAPATH,'Thres\');
    thresconfig = fullfile(thresroot,'config\');
    fn1 = fullfile(thresconfig,'settings.mat');
    fn2 = fullfile(thresconfig,'overall_failure_rate.mat');

    % Save
    save(fn1,'inpdir','resdir','allornew','preoron','thisoroverall')
    save(fn2,'overall_reject_list','overall_all_list')
catch
    warndlg([lasterr ' Forced quit.'],'Forced quit');
end

% Delete GUI
delete(h)



% ---------------------------------------------------------------------
% SETLISTSTRING subfunction - set string of listbox1
% ---------------------------------------------------------------------
function varargout = setliststring(h, eventdata, handles, varargin)

if isequal(get(handles.popupmenu2,'Value'),1)        % 'Preprocessed' mode
    sl_prepocessed(h,eventdata,handles,varargin)
elseif isequal(get(handles.popupmenu2,'Value'),2)    % 'On line' mode
    sl_online(h,eventdata,handles,varargin)
end

% ---------------------------------------------------------------------
function varargout = sl_prepocessed(h, eventdata, handles, varargin)

% Create list string from input directory
global DATAPATH
ff = handles.inpdir;
k = dir(ff);
lenlist = length(k)-2;
liststring = {};
for t = 1:lenlist
    if ~k(t+2).isdir
        liststring{end+1} = k(t+2).name(7:end-4);
    end
end

% Create list if 'New' option is set
if isequal(get(handles.popupmenu1,'Value'),2)
    k2 = dir(handles.resdir);
    lenlist2 = length(k2)-2;
    readystring = {};
    for t = 1:lenlist2
        if ~k2(t+2).isdir
            readystring{end+1} = k2(t+2).name(1:end-6);
        end
    end
    liststring = setdiff(liststring,readystring);
end

% Set list
set(handles.listbox1,'String',liststring);
guidata(h,handles)

% Reset listbox 'Value' property
if ~isempty(liststring)
    set(handles.listbox1,'Value',1)
    listbox1_Callback(h,eventdata,handles,varargin)
end

% ---------------------------------------------------------------------
function varargout = sl_online(h, eventdata, handles, varargin)



% ---------------------------------------------------------------------
% PLOTTHRES subfunction - superimpose threshold
% ---------------------------------------------------------------------
function varargout = plotthres(h, eventdata, handles, varargin)

% Set current axes
axes(handles.axes1)

% Get variables
fname = handles.fname;
segfirst = handles.segfirst;
seglast = handles.seglast;
thres = handles.threshold;
time = handles.time;

% Get axes limits
x_lim = xlim;
y_lim = ylim;

% Plot threshold
hold on
if isfield(handles,'threshandle') & ishandle(handles.threshandle)
    delete(handles.threshandle)
end
T = plot(time,thres,'g');
handles.threshandle = T;
guidata(h,handles)

% Reset axes limits
axis([x_lim(1) x_lim(2) y_lim(1) y_lim(2)]);

t = ['THRESHOLD ',fname(1:3),' ',fname(5:6),' ',num2str(segfirst),' ',num2str(seglast)];
title(t);
hold off



% ---------------------------------------------------------------------
% DISC subfunction - discrimination
% ---------------------------------------------------------------------
function vdisc = disc(unit,kuszob);

% Discriminating
disc = find(unit>=kuszob); 
discl = length(disc);
disc = [disc; unit(disc(1:discl))];
disc(1,discl+1) = length(unit) + 2;
dif = diff(disc(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j = 1:length(difn1) - 1
    [maxe,maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end