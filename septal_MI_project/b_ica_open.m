function varargout = b_ica_open(varargin)
% B_ICA_OPEN Application M-file for b_ica_open.fig
%    FIG = B_ICA_OPEN launch b_ica_open GUI.
%    B_ICA_OPEN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 01-Sep-2004 16:43:16

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    
    % Set listbox string
    global DATAPATH
    ds = [DATAPATH '\ICA\ica_newave3\list.mat'];
    load(ds)
    set(handles.listbox1,'String',list)
    
    % Store default GUI position in the handles structure for the Resize Function
    old_units = get(handles.figure1,'Units');
    set(handles.figure1,'Units','pixels');
    default_gui_pos = get(handles.figure1,'Position');
    set(handles.figure1,'Units',old_units);
    handles.gui_pos = default_gui_pos;    
    guidata(fig,handles);
    
    % Export GUI handle ('fig')
    global OPEN_GUI_HANDLE
    OPEN_GUI_HANDLE = fig;
    
    

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
% RESIZE function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h,eventdata,handles,varargin)
gui_pos = handles.gui_pos;

% Remember the old Units properties
old_units_fig = get(handles.figure1,'Units');
old_units_listbox = get(handles.listbox1,'Units');
old_units_push = get(handles.cancel,'Units');

% Set Units to pixels
set(handles.figure1,'Units','pixels');
pos = get(handles.figure1,'Position');

% Compensate if too narrow
limit = 200;
if pos(3) < limit
    if pos(1) == gui_pos(1)
        pos(3) = limit;
        set(handles.figure1,'Position',pos);
    else
        pos(1) = pos(1) + pos(3) - limit;
        pos(3) = limit;
        set(handles.figure1,'Position',pos);
    end
end

% Compensate if too low
limit = 150;
if pos(4) < limit
    if pos(2) == gui_pos(2)
        pos(4) = limit;
        set(handles.figure1,'Position',pos);
    else
        pos(2) = pos(2) + pos (4) - limit;
        pos(4) = limit;
        set(handles.figure1,'Position',pos);
    end
end

% New position of the listbox
set(handles.listbox1,'Units','pixels');
pos_list = get(handles.listbox1,'Position');
set(handles.listbox1,'Position',[pos_list(1) pos_list(2) pos(3)-gui_pos(3)+pos_list(3) pos(4)-gui_pos(4)+pos_list(4)]);
set(handles.listbox1,'Units',old_units_listbox);

% New position of the pushbuttons
set(handles.cancel,'Units','pixels');
pos_cancel = get(handles.cancel,'Position');
set(handles.cancel,'Position',[pos(3)-gui_pos(3)+pos_cancel(1) pos_cancel(2) pos_cancel(3) pos_cancel(4)]);
set(handles.cancel,'Units',old_units_push);

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% --------------------------------------------------------------------
% LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h,eventdata,handles,varargin)

% Get cell name
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cell_name = file_list{index_selected};

% Put cell name to global workspace
global CELL
CELL = cell_name;

% Get 'handles' structure
global ICA_GUI_HANDLE
h = ICA_GUI_HANDLE;
handles = guidata(h);

% Get axes handles
A = handles.axes;

% Open
for i = 1:length(handles.guidisplay)
    axes(A(i))
    if ~isempty(handles.guidisplay(i).stringinput)
        if ~isempty(handles.guidisplay(i).name)
            if ~isempty(handles.guidisplay(i).button)
                if ~isempty(handles.guidisplay(i).settings)
                    if strcmp(handles.guidisplay(i).opening_button,'current')
                        st = [''''];
                        for j = 1:length(handles.guidisplay(i).settings)
                            st = [st handles.guidisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).button) ','''...
                                handles.guidisplay(i).stringinput ''',' st ');'];
                        eval(str);
                    else
                        st = [''''];
                        for j = 1:length(handles.guidisplay(i).settings)
                            st = [st handles.guidisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).opening_button) ','''...
                                handles.guidisplay(i).stringinput ''',' st ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.guidisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                else
                    if strcmp(handles.guidisplay(i).opening_button,'current')
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).button) ','''...
                                handles.guidisplay(i).stringinput ''');'];
                        eval(str);
                    else
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).opening_button) ','''...
                                handles.guidisplay(i).stringinput ''');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.guidisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                end
            else
                if ~isempty(handles.guidisplay(i).settings)
                    st = [''''];
                    for j = 1:length(handles.guidisplay(i).settings)
                        st = [st handles.guidisplay(i).settings{j} ''','''];
                    end
                    st = st(1:end-2);
                    str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'''...
                            handles.guidisplay(i).stringinput ''',' st ');'];
                    eval(str);
                else
                    str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'''...
                            handles.guidisplay(i).stringinput ''');'];
                    eval(str);
                end
            end
        end
    else
        if ~isempty(handles.guidisplay(i).name)
            if ~isempty(handles.guidisplay(i).button)
                if ~isempty(handles.guidisplay(i).settings)
                    if strcmp(handles.guidisplay(i).opening_button,'current')
                        st = [''''];
                        for j = 1:length(handles.guidisplay(i).settings)
                            st = [st handles.guidisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).button) ',' st ');'];
                        eval(str);
                    else
                        st = [''''];
                        for j = 1:length(handles.guidisplay(i).settings)
                            st = [st handles.guidisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).opening_button) ',' st ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.guidisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                else
                    if strcmp(handles.guidisplay(i).opening_button,'current')
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).button) ');'];
                        eval(str);
                    else
                        str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,'...
                                num2str(handles.guidisplay(i).opening_button) ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.guidisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                end
            else
                if ~isempty(handles.guidisplay(i).settings)
                    st = [''''];
                    for j = 1:length(handles.guidisplay(i).settings)
                        st = [st handles.guidisplay(i).settings{j} ''','''];
                    end
                    st = st(1:end-2);
                    str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles,' st ');'];
                    eval(str);
                else
                    str = ['b_ica_gui4(''' handles.guidisplay(i).name ''',h,handles);'];
                    eval(str);
                end
            end
        end 
    end
end

for i = 1:length(handles.figdisplay)
    figure(handles.figure(i))
%     delete(get(handles.figure(i),'Children'))
    if ~isempty(handles.figdisplay(i).stringinput)
        if ~isempty(handles.figdisplay(i).name)
            if ~isempty(handles.figdisplay(i).button)
                if ~isempty(handles.figdisplay(i).settings)
                    if strcmp(handles.figdisplay(i).opening_button,'current')
                        st = [''''];
                        for j = 1:length(handles.figdisplay(i).settings)
                            st = [st handles.figdisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).button) ','...
                                handles.figdisplay(i).stringinput ',' st ');'];
                        eval(str);
                    else
                        st = [''''];
                        for j = 1:length(handles.figdisplay(i).settings)
                            st = [st handles.figdisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).opening_button) ','...
                                handles.figdisplay(i).stringinput ',' st ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.figdisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                else
                    if strcmp(handles.figdisplay(i).opening_button,'current')
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).button) ','...
                                handles.figdisplay(i).stringinput ');'];
                        eval(str);
                    else
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).opening_button) ','...
                                handles.figdisplay(i).stringinput ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.figdisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                end
            else
                if ~isempty(handles.figdisplay(i).settings)
                    st = [''''];
                    for j = 1:length(handles.figdisplay(i).settings)
                        st = [st handles.figdisplay(i).settings{j} ''','''];
                    end
                    st = st(1:end-2);
                    str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                            handles.figdisplay(i).stringinput ',' st ');'];
                    eval(str);
                else
                    str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                            handles.figdisplay(i).stringinput ');'];
                    eval(str);
                end
            end
        end
    else
        if ~isempty(handles.figdisplay(i).name)
            if ~isempty(handles.figdisplay(i).button)
                if ~isempty(handles.figdisplay(i).settings)
                    if strcmp(handles.figdisplay(i).opening_button,'current')
                        st = [''''];
                        for j = 1:length(handles.figdisplay(i).settings)
                            st = [st handles.figdisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).button) ',' st ');'];
                        eval(str);
                    else
                        st = [''''];
                        for j = 1:length(handles.figdisplay(i).settings)
                            st = [st handles.figdisplay(i).settings{j} ''','''];
                        end
                        st = st(1:end-2);
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).opening_button) ',' st ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.figdisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                else
                    if strcmp(handles.figdisplay(i).opening_button,'current')
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).button) ');'];
                        eval(str);
                    else
                        str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,'...
                                num2str(handles.figdisplay(i).opening_button) ');'];
                        eval(str);
                        B = b_ica_gui4('get_button_handle',h,handles,gca,handles.figdisplay(i).opening_button);
                        b_ica_gui4('set_font_weight',h,handles,gca,B)
                    end
                end
            else
                if ~isempty(handles.figdisplay(i).settings)
                    st = [''''];
                    for j = 1:length(handles.figdisplay(i).settings)
                        st = [st handles.figdisplay(i).settings{j} ''','''];
                    end
                    st = st(1:end-2);
                    str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles,' st ');'];
                    eval(str);
                else
                    str = ['b_ica_gui4(''' handles.figdisplay(i).name ''',h,handles);'];
                    eval(str);
                end
            end
        end 
    end
end

% Let the GUI be the current figure
figure(handles.figure1)



% --------------------------------------------------------------------
% CANCEL button callback - delete GUI
% --------------------------------------------------------------------
function varargout = cancel_Callback(h,eventdata,handles,varargin)

% Delete GUI
delete(handles.figure1)