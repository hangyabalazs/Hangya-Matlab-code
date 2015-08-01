function varargout = b_hcn_gui(varargin)
%HCN_GUI    Displays analysed segments and burst parameters.
%   HCN_GUI shows analysed segments with bursts crosslined. Intraburst spikes
%   are red, extraburst spikes are blue. The minimum of inter-first-spike
%   interval CV (coefficient of variation) values from 2 to 7 clusters, burst
%   length CV and intraburst interval CV are displayed in the edit text box.
%
%   See also HCN_ANALYSIS and ITCLUST2.

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
	end

end




% --------------------------------------------------------------------
% LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)

% Load
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');

fln = file_list{index_selected};
fs = findstr(fln,'_');
fname = fln(fs(2)+1:fs(4)-1);
segfirst = str2num(fln(fs(4)+1:fs(5)-1));
seglast = str2num(fln(fs(5)+1:end-4));

ff2 = [DATAPATH,'HCN\Analysis2\HCN_ICA_',fname,'_',num2str(segfirst),'_',num2str(seglast)];   % load variance vectors
load(ff2)

ff = fullfile(DATAPATH,'HCN\Analysis2',file_list{index_selected});   % load 'Burst' matrix
load(ff)

% where2 = ['f:\raw_data\hcn\all\'];;
% files2 = dir(where2);
% ns_long = cell(1,length(files2));
% ns_short = cell(1,length(files2));
% for i = 3:length(files2)
%     ns_long{i} = files2(i).name;
%     ns_short{i} = ns_long{i}(1:6);
% end
% scp = strcmp(fname,ns_short);
% fnd = find(scp);
% fln = ns_long{fnd};
% ff2 = [where2 fln];
% data = b_load_data(ff2);   % load data
% unit = data(segfirst:seglast,2);
% unit = unit';
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
axes(handles.axes1)
bdf1 = get(handles.axes1,'ButtonDownFcn');
z = zeros(1,length(time));
z(vdisc) = 1;
z(vdisc+1) = -1;
plot1_handle = plot(time,z,'b');
axis([time(1) time(end) -1.5 1.5])
hold on;
sbd = size(Burst,2);
hans = zeros(1,sbd);
for j = 1:sbd
    ind1 = vdisc(Burst(1,j));
    ind2 = vdisc(Burst(2,j));
    hans(j) = plot(time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r');
    rajz = [time(vdisc(Burst(1,j))) time(vdisc(Burst(2,j))); 0.2 0.7];
    line_handle = line(rajz(1,:),rajz(2,:),'Color','k');
end
hold off

% Set edit string
mfsp = min(FirstSpikeCv(2:7));
im = find(FirstSpikeCv==mfsp);
bl = BurstLengthCv(im(1));
ib = IntraBurstIvCv(im(1));
mfsp_bl_ib = [num2str(mfsp) ';   ' num2str(bl) ';   ' num2str(ib)];
set(handles.edit1,'String',mfsp_bl_ib)

% Plot on axis1 - unit with bursts
% axes(handles.axes1)
% bdf1 = get(handles.axes1,'ButtonDownFcn');
% plot3(mfsp,bl,ib);

% Set buttondown functions
% handles.plot2_handle = plot2_handle;
% set(handles.plot2_handle,'ButtonDownFcn',bdf2)
% set(handles.axes2,'ButtonDownFcn',bdf2);

handles.plot1_handle = plot1_handle;
set(handles.plot1_handle,'ButtonDownFcn',bdf1)
set(handles.axes1,'ButtonDownFcn',bdf1);



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
    'Position',[pos(3)-gui_pos(3)+pos_list(1) pos(4)-gui_pos(4)+pos_list(2) pos_list(3) pos_list(4)]);

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
        listbox1_Callback(h, eventdata, handles, varargin);
    end
otherwise
    if index_selected - 1 >= 1
        set(handles.listbox1,'Value',index_selected-1)
        listbox1_Callback(h, eventdata, handles, varargin);
    end
end



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
        lb2(h, eventdata, handles, varargin);
    end
otherwise
    if index_selected - 1 >= 1
        set(handles.listbox1,'Value',index_selected-1)
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
ff = fullfile(DATAPATH,'HCN\Analysis2');
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


% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)

