function varargout = b_2entropy_gui(varargin)
%

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'new');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
    
    % Input directory, wavelet scale vector
    global DATAPATH
    handles.inpdir = get(handles.edit1,'String');
    load([handles.inpdir 'settings\ScaleVector']);
    handles.scalevector = f;
    guidata(fig,handles)
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
    set(handles.checkbox6,'Value',1)
    set(handles.checkbox7,'Value',1)
    set(handles.checkbox8,'Value',1)
    set(handles.checkbox9,'Value',1)
    

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
set(handles.listbox1,'Units','pixels');
pos_list = get(handles.listbox1,'Position');
set(handles.listbox1,...
    'Position',[pos(3)-gui_pos(3)+pos_list(1) pos_list(2) pos_list(3) pos(4)-gui_pos(4)+pos_list(4)]);
set(handles.listbox1,'Units',old_units_listbox);

% New position of the edit text
pos_edit = cell(1,1);
set(handles.edit1,'Units','pixels');
pos_edit{1} = get(handles.edit1,'Position');
set(handles.edit1,...
    'Position',[pos_edit{1}(1) pos_edit{1}(2) pos_edit{1}(3) pos_edit{1}(4)]);
set(handles.edit1,'Units',old_units_edit);

% New position of the static text above the listbox
pos_text = cell(1,4);
set(handles.text1,'Units','pixels');
pos_text{1} = get(handles.text1,'Position');
set(handles.text1,...
    'Position',[pos(3)-gui_pos(3)+pos_text{1}(1) pos(4)-gui_pos(4)+pos_text{1}(2) pos_text{1}(3) pos_text{1}(4)]);
set(handles.text1,'Units',old_units_text);

% New position of the static text in the 'input directory' panel
set(handles.text2,'Units','pixels');
pos_text{2} = get(handles.text2,'Position');
set(handles.text2,...
    'Position',[pos_text{2}(1) pos_text{2}(2) pos_text{2}(3) pos_text{2}(4)]);
set(handles.text2,'Units',old_units_text);

% New position of the axes
pos_axes = cell(1,2);
set(handles.axes1,'Units','pixels');
pos_axes{1} = get(handles.axes1,'Position');
set(handles.axes1,...
    'Position',[pos_axes{1}(1) pos_axes{1}(2)+((pos(4)-gui_pos(4))/2)...
        pos(3)-gui_pos(3)+pos_axes{1}(3) ((pos(4)-gui_pos(4))/2)+pos_axes{1}(4)]);
set(handles.axes1,'Units',old_units_axes);
set(handles.axes2,'Units','pixels');
pos_axes{2} = get(handles.axes2,'Position');
set(handles.axes2,...
    'Position',[pos_axes{2}(1) pos_axes{2}(2)...
        pos(3)-gui_pos(3)+pos_axes{2}(3) ((pos(4)-gui_pos(4))/2)+pos_axes{2}(4)]);
set(handles.axes2,'Units',old_units_axes);

% New position of the checkbox
pos_check = cell(1,9);
for t = 1:9
    eval(['set(handles.checkbox',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_check{t} = get(handles.checkbox',int2str(t),',''Position'');']);
    eval(['set(handles.checkbox',int2str(t),...
            ',''Position'',[pos_check{t}(1) pos_check{t}(2) pos_check{t}(3) pos_check{t}(4)]);']);
    eval(['set(handles.checkbox',int2str(t),',''Units'',old_units_check)']);
end

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% --------------------------------------------------------------------
% LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)

% Load
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
load(fullfile(handles.inpdir,file_list{index_selected}))

% Time
time = [30:0.5:120];
time = time(1:length(Rhyabs));

% Plot on axis1 - uncertainity coefficients
axes(handles.axes1)
plot(time,Ruxyabs,'b.-','markersize',12)
hold on
plot(time,Ruyxabs,'r.-','markersize',12)
hold off
legend('U(eeg->unit)','U(unit->eeg)');

% Plot on axis2 - entropy values
axes(handles.axes2)
is1 = get(handles.checkbox1,'Value');
is2 = get(handles.checkbox2,'Value');
is3 = get(handles.checkbox3,'Value');
is4 = get(handles.checkbox4,'Value');
is5 = get(handles.checkbox5,'Value');
is6 = get(handles.checkbox6,'Value');
is7 = get(handles.checkbox6,'Value');
is8 = get(handles.checkbox6,'Value');
is9 = get(handles.checkbox6,'Value');
legend_matrix = {};
if is1
    plot(time,Rhyabs,'b.-','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(eeg)';
end
if is2
    plot(time,Rhxabs,'r.-','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(unit)';
end
if is3
    plot(time,Rhxyabs,'g.-','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(eeg,unit)';
end
if is4
    plot(time,Rixyabs,'g.:','markersize',12)
    hold on
    legend_matrix{end+1} = 'I(eeg&unit)';
end
if is5
    plot(time,Rixynormabs,'g.:','markersize',12)
    hold on
    legend_matrix{end+1} = 'I(eeg&unit) normalized';
end
if is6
    plot(time,Rhxcyabs,'r.--','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(unit|eeg)';
end
if is7
    plot(time,Rhycxabs,'b.--','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(eeg|unit)';
end
if is8
    plot(time,Rrelshanyabs,'b.-.','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(eeg) relative';
end
if is9
    plot(time,Rrelshanxabs,'r.-.','markersize',12)
    hold on
    legend_matrix{end+1} = 'H(unit) relative';
end
hold off
legend(legend_matrix);

% Set buttondown functions
bdf = 'b_2entropy_gui(''axes_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
ch1 = get(handles.axes1,'Children');
set(ch1,'ButtonDownFcn',bdf)
set(handles.axes1,'ButtonDownFcn',bdf);
ch2 = get(handles.axes2,'Children');
set(ch2,'ButtonDownFcn',bdf)
set(handles.axes2,'ButtonDownFcn',bdf);



% --------------------------------------------------------------------
% POPUPMENU1 callback - switch mode
% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)

% Enable
pop = get(handles.popupmenu1,'Value');
if isequal(pop,1)
    set(handles.popupmenu2,'Enable','off')
    set(handles.popupmenu3,'Enable','off')
    set(handles.checkbox1,'Enable','on')
    set(handles.checkbox2,'Enable','on')
    set(handles.checkbox3,'Enable','on')
    set(handles.checkbox4,'Enable','on')
    set(handles.checkbox5,'Enable','on')
    set(handles.checkbox6,'Enable','on')
    set(handles.checkbox7,'Enable','on')
    set(handles.checkbox8,'Enable','on')
    set(handles.checkbox9,'Enable','on')
elseif isequal(pop,2)
    set(handles.popupmenu2,'Enable','on')
    set(handles.popupmenu3,'Enable','on')
    set(handles.checkbox1,'Enable','off')
    set(handles.checkbox2,'Enable','off')
    set(handles.checkbox3,'Enable','off')
    set(handles.checkbox4,'Enable','off')
    set(handles.checkbox5,'Enable','off')
    set(handles.checkbox6,'Enable','off')
    set(handles.checkbox7,'Enable','off')
    set(handles.checkbox8,'Enable','off')
    set(handles.checkbox9,'Enable','off')
end

% Refresh
listbox1_Callback(h, eventdata, handles, varargin);



% --------------------------------------------------------------------
% POPUPMENU2 callback - switch image on axes1
% --------------------------------------------------------------------
function varargout = popupmenu2_Callback(h, eventdata, handles, varargin)

% Load
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
load(fullfile(handles.inpdir,file_list{index_selected}))

% Plot on axis1
axes(handles.axes1)
switch get(handles.popupmenu2,'Value')
case 1
    imagesc(Ruxyabs)
case 2
    imagesc(Ruyxabs)
case 3
    imagesc(Rhyabs)
case 4
    imagesc(Rhxabs)
case 5
    imagesc(Rhxyabs)
case 6
    imagesc(Rixyabs)
case 7
    imagesc(Rixynormabs)
case 8
    imagesc(Rhxcyabs)
case 9
    imagesc(Rhycxabs)
case 10
    imagesc(Rrelshanabs)
case 11
    msgbox('In future versions only.')
end



% --------------------------------------------------------------------
% POPUPMENU3 callback - switch image on axes2
% --------------------------------------------------------------------
function varargout = popupmenu3_Callback(h, eventdata, handles, varargin)

% Load
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
load(fullfile(handles.inpdir,file_list{index_selected}))

% Plot on axis1
axes(handles.axes2)
switch get(handles.popupmenu3,'Value')
case 1
    imagesc(Ruxyabs)
case 2
    imagesc(Ruyxabs)
case 3
    imagesc(Rhyabs)
case 4
    imagesc(Rhxabs)
case 5
    imagesc(Rhxyabs)
case 6
    imagesc(Rixyabs)
case 7
    imagesc(Rixynormabs)
case 8
    imagesc(Rhxcyabs)
case 9
    imagesc(Rhycxabs)
case 10
    imagesc(Rrelshanabs)
case 11
    msgbox('In future versions only.')
end



% --------------------------------------------------------------------
% EDIT TEXT1 callback - input directory
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
inpdir = get(handles.edit1,'String');
if ~isdir(inpdir)
    errordlg('Input directory do not exist.')
else
    handles.inpdir = inpdir;
end
guidata(h,handles);
setliststring(h,eventdata,handles,varargin);



% --------------------------------------------------------------------
% AXES BUTTONDOWN function - step one file
% --------------------------------------------------------------------
function varargout = axes_ButtonDownFcn(h, eventdata, handles, varargin)
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
ff = handles.inpdir;
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
function varargout = checkbox7_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox8_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox9_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)