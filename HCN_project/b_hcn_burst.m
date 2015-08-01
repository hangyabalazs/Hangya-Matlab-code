function varargout = b_hcn_burst(varargin)
%HCN_BURST  Semiautomatic burst selector for HCN project.
%   HCN_BURST loads iterative cluster analysis results (see ITCLUST for details).
%   User can set the number of clusters from 2 to 20 in the listbox. Variability
%   analysis results (variance plot) is getting displayed in the upper axes, while
%   the lower axes shows the unit with bursts.
%   Optimal cluster number can be set by editing the edit box and pressing enter.
%   After this action, HCN_BURST GUI steps to the next registration segment in the
%   input directory (which user can modify only through editing the program code).
%   For nonbursty cells, you have to assign 0 as optimal cluster number.
%
%   Optimal clusternumber is getting saved for HCN_ANALYSIS. Bursty and nonbursty
%   segments' input files are sorted into different subdirectories.
%
%   See also ITCLUST and HCN_ANALYSIS.


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
    
    % Maximal cluster number
    handles.max_clusno = 20;
    
    % Assign first cell as 'fln' (filename)
    global DATAPATH
    files = dir([DATAPATH 'HCN\Burst']);
    files = files(3:end);
    files2 = {};
    for i = 1:length(files)
        if ~files(i).isdir
            files2{end+1} = files(i).name;
        end
    end
    handles.fln = files2{1};
    handles.flnno = 1;
    handles.files = files2;
    
    % Set string of listbox
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

% Load ICA data
global DATAPATH
fln = handles.fln;
fs = findstr(fln,'_');
fname = fln(fs(2)+1:fs(4)-1);
segfirst = str2num(fln(fs(4)+1:fs(5)-1));
seglast = str2num(fln(fs(5)+1:end-4));
ff = [DATAPATH 'HCN\Burst\' fln];
load(ff)

% Plot on axis1 - variance plot
axes(handles.axes1)
d = length(IntraBurstIvCv);

isintra = get(handles.checkbox1,'Value');
isextra = get(handles.checkbox2,'Value');
isinter = get(handles.checkbox3,'Value');
isfirst = get(handles.checkbox4,'Value');
isallfirst = get(handles.checkbox5,'Value');

if isintra
    plot1_handle = plot([1:d],IntraBurstIvCv,'g');
    hold on
end
if isextra
    plot1_handle = plot([1:d],ExtraBurstIvCv,'b');
    hold on
end
if isinter
    plot1_handle = plot([1:d],InterBurstIvCv,'c');
    hold on
end
if isfirst
    plot1_handle = plot([1:d],FirstSpikeCv,'m');
    hold on
end
if isallfirst
    plot1_handle = plot([1:d],AllFirstSpikeCv,'k');
    hold on
end

legend('intraburstivcv','extraburstivcv','interburstivcv','firstspikecv','allfirstspikecv',2);
xlim([1 d]);
x_lim = xlim;
y_lim = ylim;
kk = (y_lim(2) - y_lim(1)) / 3;
length_vdisc = length(Vdisc);
if length_vdisc < 100,
    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length_vdisc)],'Color','r');
else
    text(x_lim(1)+1,y_lim(2)-kk,['{\itNumber of spikes : }',num2str(length_vdisc)],'Color','k');
end    
hold off

% Create appropriate 'vdisc'
if segfirst > 2000000
    pseudoindex = segfirst - 2000000;
    vdisc = Vdisc - pseudoindex;
else
    vdisc = Vdisc - segfirst;
end

% Create 'time' vector
dt = 0.0001;
len = seglast - segfirst;
time = [1:len] * dt;

% Plot on axis1 - unit with bursts
index_selected = get(handles.listbox1,'Value');
dec = index_selected + 1;
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
z = zeros(1,length(time));
z(vdisc) = 1;
z(vdisc+1) = -1;
plot1_handle = plot(time,z,'b');
axis([time(1) time(end) -1.5 1.5])
hold on;
sbd = size(Burst{dec},2);
hans = zeros(1,sbd);
for j = 1:sbd
    ind1 = vdisc(Burst{dec}(1,j));
    ind2 = vdisc(Burst{dec}(2,j));
    hans(j) = plot(time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r');
    rajz = [time(vdisc(Burst{dec}(1,j))) time(vdisc(Burst{dec}(2,j))); 0.2 0.7];
    line_handle = line(rajz(1,:),rajz(2,:),'Color','k');
end
hold off

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
set(handles.axes2,'ButtonDownFcn',bdf2);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

% Refresh 'handles' structure
guidata(h,handles);



% --------------------------------------------------------------------
% RESIZE function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
gui_pos = handles.gui_pos;

% Remember the old Units properties
old_units_fig = get(handles.figure1,'Units');
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
set(handles.listbox1,'Units','pixels');
pos_list = get(handles.listbox1,'Position');
set(handles.listbox1,...
    'Position',[pos(3)-gui_pos(3)+pos_list(1) pos_list(2) pos_list(3) pos(4)-gui_pos(4)+pos_list(4)]);

% New position of the checkbox
pos_check = cell(1,5);
for t = 1:5
    eval(['set(handles.checkbox',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_check{t} = get(handles.checkbox',int2str(t),',''Position'');']);
    eval(['set(handles.checkbox',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_check{t}(1) pos_check{t}(2) pos_check{t}(3) pos_check{t}(4)]);']);
end

% New position of the edit text
set(handles.edit1,'Units','pixels');
pos_edit = get(handles.edit1,'Position');
set(handles.edit1,...
    'Position',[pos_edit(1) pos_edit(2) pos(3)-gui_pos(3)+pos_edit(3) pos_edit(4)]);

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

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);
set(handles.axes1,'Units',old_units_axes);
set(handles.listbox1,'Units',old_units_listbox);
set(handles.edit1,'Units',old_units_edit);
set(handles.checkbox1,'Units',old_units_check);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% --------------------------------------------------------------------
% EDIT TEXT CALLBACK - switch to next cell
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)

% Get 'string'
opt_clusno = get(handles.edit1,'String');       % get optimal cluster number
opt_clusno = str2num(opt_clusno);
fln = handles.fln;
fs = findstr(fln,'_');
fname = fln(fs(2)+1:fs(4)-1);
segfirst = str2num(fln(fs(4)+1:fs(5)-1));
seglast = str2num(fln(fs(5)+1:end-4));
newfn = ['OPTIMAL_CLUSNO_' fname '_' num2str(segfirst) '_' num2str(seglast) '.mat'];

% Change directory
mmm = pwd;
global DATAPATH
cd([DATAPATH,'HCN\Burst\']); 

% Save
respath1 = [DATAPATH,'HCN\Burst\bursty\'];
respath2 = [DATAPATH,'HCN\Burst\nonbursty'];
newpth = [respath1 newfn];
if isequal(opt_clusno,0)
    eval(['copyfile(''',handles.fln,''',''',respath2,''');']);
else
    eval(['copyfile(''',handles.fln,''',''',respath1,''');']);
    eval(['save ' newpth ' opt_clusno']);
end

% Increase 'flnno' (filename number) by one
handles.flnno = handles.flnno + 1;

% New 'fln' (filename)
handles.fln = handles.files{handles.flnno};

% Refresh 'handles' structure
guidata(h,handles)

% Plot new file
lb(h,eventdata,handles,varargin)

% Reset directory
cd(mmm)



% --------------------------------------------------------------------
% AXES2 BUTTONDOWN function - step one figure
% --------------------------------------------------------------------
function varargout = axes2_ButtonDownFcn(h, eventdata, handles, varargin)
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



% ---------------------------------------------------------------------
% SETLISTSTRING subfunction - set string of listbox
% ---------------------------------------------------------------------
function varargout = setliststring(h, eventdata, handles, varargin)

% Get number of clusters
d = handles.max_clusno;

% Create list string
list = cell(1,d-1);
for t = 2:d
    list{t-1} = ['Number of clusters = ',int2str(t)];
end

% Set list string
set(handles.listbox1,'String',list);