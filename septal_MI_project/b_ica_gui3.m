function varargout = b_ica_gui3(varargin)
%ICA_GUI3   Graphical User Interface for Iterative Cluster Analysis.
%   The only difference between ICA_GUI2B and ICA_GUI3 is that in ICA_GUI3
%   you can switch between iterative and natural clustering method in
%   'Options' - 'Cloustering' menu.
%
%   Note: not debugged!
%
%   See also ICA_GUI2B and NATCLUST.

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
    
    % Set default value of whichplot
    handles.whichplot = 'variance';
    
    % Set default value of whichperiodogram
    global WHICHPERIODOGRAM
    WHICHPERIODOGRAM = 1;
    
    % Set default value of varcell
    handles.varcell = [];
    
    % Set default value of accessoric figure handles
    handles.accessoric_axis_handle= [];
    handles.accessoric_fig_handle= [];
    handles.subplot1_handle = [];
    handles.subplot2_handle = [];
    handles.greyplot1_handle = [];
    handles.greyplot2_handle = [];
    handles.greyplot3_handle = [];
    handles.greyplot4_handle = [];
    handles.line_handles = [];
    
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
% ALLVAR LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles, varargin)
switch handles.whichplot
case 'variance'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb(h, eventdata, handles, varargin)
    end
case 'waveletvector'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb_wavevec(h, eventdata, handles, varargin)
    end
case 'lombperiodogram'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb_lombper(h, eventdata, handles, varargin)
    end
case 'fft'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb_fft(h, eventdata, handles, varargin)
    end
end

% --------------------------------------------------------------------
function varargout = lb(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
if isequal(handles.setinput,['base'])
    nm = '';
else
    nm = ['ica_beta_' handles.setinput];
end
cl = file_list{index_selected}(11:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\',nm];
fln = ['THETA_ICA_' cl];
ff = fullfile(pth,fln);
load(ff)
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

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

plot1_handle = [];
if isintra
    plot1_handle(end+1) = plot([1:d],IntraBurstIvVar,'g');
    hold on
end
if isextra
    plot1_handle(end+1) = plot([1:d],ExtraBurstIvVar,'b');
    hold on
end
if isinter
    plot1_handle(end+1) = plot([1:d],InterBurstIvVar,'c');
    hold on
end
if isfirst
    plot1_handle(end+1) = plot([1:d],FirstSpikeVar,'m');
    hold on
end
if isallfirst
    plot1_handle(end+1) = plot([1:d],AllFirstSpikeVar,'k');
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
end        

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
index_selected = get(handles.listbox2,'Value');
if isequal(index_selected,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
    hans = [];
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected + 1;
    plot2_handle = plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    sbd = size(Burst{dec},2);
    hans = zeros(1,sbd);
    for j = 1:sbd
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        hans(j) = plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r');
%         rajz = [Time(Vdisc_new(Burst{dec}(1,j))) Time(Vdisc_new(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
handles.plot1_handle = plot1_handle;
set(handles.axes1,'ButtonDownFcn',bdf1);

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
if ~ishandle(h)
    h = handles.figure1;
end
guidata(h,handles);

% Set keypress function
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

% --------------------------------------------------------------------
function lb_wavevec(h, eventdata, handles, varargin)

% Load wavelet vectors
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave\',file_list{index_selected});
load(ff)

% Load ICA data
cl = file_list{index_selected}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
ff = fullfile(pth,fln);
load(ff)
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

% Plot on axis1 - wavelet averege vectors
axes(handles.axes1)
bdf1 = get(handles.axes1,'ButtonDownFcn');
plot1_handle = plot(Wavevec(1,:),'r');
hold on
plot1_handle = plot(Wavevec(2,:),'b');
plot1_handle = plot(Wavevec(3,:),'g');
plot1_handle = plot(Wavevec(4,:),'m');
legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
hold off
y_lim = ylim;
intlen =  length(Wavevec);
axis([0 intlen y_lim(1) y_lim(2)]);

% Create list for listbox2
d = length(IntraBurstIvVar);
list = cell(1,d-1);
list{1} = 'Return Plot';
for t = 2:d-1
    list{t} = ['Number of clusters +1 = ',int2str(t+1)];
end
set(handles.listbox2,'String',list);

% Set edit string
set(handles.edit1,'String','WAVELET VECTORS')

% Plot on axis2 - return plot or return plot or spike train with structural elements
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
index_selected = get(handles.listbox2,'Value');
if isequal(index_selected,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected + 1;
    plot2_handle = plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
handles.plot1_handle = plot1_handle;
set(handles.axes1,'ButtonDownFcn',bdf1);

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
if ~ishandle(h)
    h = handles.figure1;
end
guidata(h,handles)

% Set keypress function
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

% --------------------------------------------------------------------
function lb_lombper(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cl = file_list{index_selected}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
ff = fullfile(pth,fln);
load(ff)
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

d = length(IntraBurstIvVar);
MaxClusterNumber = d;
Vdisc_new = Vdisc - (i_first - 320000);

% Plot on axis1 - wavelet averege vectors or Lomb periodogram
axes(handles.axes1)
index_selected = get(handles.listbox2,'Value');
bdf1 = get(handles.axes1,'ButtonDownFcn');
if isequal(index_selected,1)
    index_selected_lb1 = get(handles.listbox1,'Value');
    file_list = get(handles.listbox1,'String');
    ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave\',file_list{index_selected_lb1});
    load(ff)
    plot1_handle = plot(Wavevec(1,:),'r');
    hold on
    plot1_handle = plot(Wavevec(2,:),'b');
    plot1_handle = plot(Wavevec(3,:),'g');
    plot1_handle = plot(Wavevec(4,:),'m');
    legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
    hold off
    y_lim = ylim;
    intlen =  length(Wavevec);
    axis([0 intlen y_lim(1) y_lim(2)]);
else
    dec = index_selected + 1;
    bas = [];
    for bno = 1:size(Burst{dec},2)
        b = Vdisc_new(Burst{dec}(1,bno):Burst{dec}(2,bno));
        bas = [bas b];  %bas is for "burst all spike" and contains the locations of all intraburst spikes
    end
%     fbs = Vdisc_new(Burst{dec}(1,:));
%     lenu = length(Time);
%     vdisc_for_lomb = Vdisc_new;
    fname = cl(1:6);
    period_handle = b_lombper_for_ica(fname,i_first,i_second,dec);  %Lomb periodogram is calculated for the
    plot1_handle = period_handle;                                   %"burst all spike" array
%     [pxx,pyy,jmax,prob,z,effm,period_handle] = b_lomb_period(vdisc_for_lomb,bas,lenu);
end

% Create list for listbox2
list = cell(1,d-1);
list{1} = 'Return Plot';
for t = 2:d-1
    list{t} = ['Number of clusters +1 = ',int2str(t+1)];
end
set(handles.listbox2,'String',list);

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
index_selected = get(handles.listbox2,'Value');
if isequal(index_selected,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected + 1;
    plot2_handle = plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%         rajz = [Time(Vdisc_new(Burst{dec}(1,j))) Time(Vdisc_new(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
handles.plot1_handle = plot1_handle;
set(handles.axes1,'ButtonDownFcn',bdf1);

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
if ~ishandle(h)
    h = handles.figure1;
end
guidata(h,handles);

% Set keypress function
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

% --------------------------------------------------------------------
function varargout = lb_fft(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected1 = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cl = file_list{index_selected1}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
ff = fullfile(pth,fln);
load(ff)
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

d = length(IntraBurstIvVar);
MaxClusterNumber = d;

% Plot on axis1 - wavelet magnitude average vectors' Fourier transform
axes(handles.axes1)
bdf1 = get(handles.axes1,'ButtonDownFcn');
index_selected2 = get(handles.listbox2,'Value');
file_list = get(handles.listbox1,'String');
ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave\',file_list{index_selected1});
load(ff)
lenwv = size(Wavevec,2);
pp1 = log2(lenwv);
pp2 = floor(pp1);
pp3 = 2^pp2;
global WHICHPERIODOGRAM
whichperiod = WHICHPERIODOGRAM;
switch whichperiod
case 1
%     [p_1 f_1] = periodogram(Wavevec(4,:),[],pp3,100);
%     k_1 = f_1 * 100 / pp3;
%     plot1_handle = plot(k_1(4:end),p_1(4:end),'m');
%     ym = max([p_1(5:end)]);
%     axis([0 5 0 2*ym]);
%     text(3,1.5*ym,'{\it1 - 3 Hz}')
    [p_1 f_1] = pwelch(Wavevec(4,:),[],[],pp3,100);
    p_11 = [0;p_1(8:end);p_1(end)];
    f_11 = [0;f_1(8:end);0];
    plot1_handle = fill(f_11,p_11,'m');
    ym = max([p_1(8:end)]);
    text(25,0.5*ym,'{\it1 - 3 Hz}')
case 2
%     [p_2 f_2] = periodogram(Wavevec(3,:),[],pp3,100);
%     k_2 = f_2 * 100 / pp3;
%     plot1_handle = plot(k_2(4:end),p_2(4:end),'g');
%     ym = max([p_2(5:end)]);
%     axis([0 5 0 2*ym]);
%     text(3,1.5*ym,'{\it3 - 6 Hz}')
    [p_2 f_2] = pwelch(Wavevec(3,:),[],[],pp3,100);
    p_22 = [0;p_2(8:end);p_2(end)];
    f_22 = [0;f_2(8:end);0];
    plot1_handle = fill(f_22,p_22,'g');
    ym = max([p_2(8:end)]);
    text(25,0.5*ym,'{\it3 - 6 Hz}')
case 3
%     [p_3 f_3] = periodogram(Wavevec(2,:),[],pp3,100);
%     k_3 = f_3 * 100 / pp3;
%     plot1_handle = plot(k_3(4:end),p_3(4:end),'b');
%     ym = max([p_3(5:end)]);
%     axis([0 5 0 2*ym]);
%     text(3,1.5*ym,'{\it6 - 20 Hz}')
    [p_3 f_3] = pwelch(Wavevec(2,:),[],[],pp3,100);
    p_33 = [0;p_3(8:end);p_3(end)];
    f_33 = [0;f_3(8:end);0];
    plot1_handle = fill(f_33,p_33,'b');
    ym = max([p_3(8:end)]);
    text(25,0.5*ym,'{\it6 - 20 Hz}')
case 4
%     [p_4 f_4] = periodogram(Wavevec(1,:),[],pp3,100);
%     k_4 = f_4 * 100 / pp3;
%     plot1_handle = plot(k_4(4:end),p_4(4:end),'r');
%     ym = max([p_4(5:end)]);
%     axis([0 5 0 2*ym]);
%     text(3,1.5*ym,'{\it20 - 50 Hz}')
    [p_4 f_4] = pwelch(Wavevec(1,:),[],[],pp3,100);
    p_44 = [0;p_4(8:end);p_4(end)];
    f_44 = [0;f_4(8:end);0];
    plot1_handle = fill(f_44,p_44,'r');
    ym = max([p_4(8:end)]);
    text(25,0.5*ym,'{\it20 - 50 Hz}')
end

% Create list for listbox2
list = cell(1,d-1);
list{1} = 'Return Plot';
for t = 2:d-1
    list{t} = ['Number of clusters +1 = ',int2str(t+1)];
end
set(handles.listbox2,'String',list);

% Plot on axis2 - return plot or Lomb periodogram
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
if isequal(index_selected2,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    Vdisc_new = Vdisc - (i_first - 320000);
    dec = index_selected2 + 1;
    bas = [];
    for bno = 1:size(Burst{dec},2)
        b = Vdisc_new(Burst{dec}(1,bno):Burst{dec}(2,bno));
        bas = [bas b];
    end
%     fbs = Vdisc_new(Burst{dec}(1,:));
%     lenu = length(Time);
%     vdisc_for_lomb = Vdisc_new;
    fname = cl(1:6);
    period_handle = b_lombper_for_ica(fname,i_first,i_second,dec);
    plot2_handle = period_handle;
%     [pxx,pyy,jmax,prob,z,effm,period_handle] = b_lomb_period(vdisc_for_lomb,bas,lenu);
end

% Plot on the accessoric figure - Wavelet average vectors on subplot1
handles.greyplot1_handle = [];
handles.greyplot2_handle = [];
handles.greyplot3_handle = [];
handles.greyplot4_handle = [];
handles.line_handles = [];
if ~isempty(handles.accessoric_axis_handle)
    H = handles.accessoric_axis_handle;
    F = handles.accessoric_fig_handle;
    S1 = handles.subplot1_handle;
    S2 = handles.subplot2_handle;
    axes(H(1));
    plot(Wavevec(1,:),'r');
    hold on
    plot(Wavevec(2,:),'b');
    plot(Wavevec(3,:),'g');
    plot(Wavevec(4,:),'m');
    legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
    hold off
    y_lim = ylim;
    intlen =  length(Wavevec);
    axis([0 intlen y_lim(1) y_lim(2)]);
    H(1) = get(gcf,'CurrentAxes');
    axes(H(2));
else
    F = figure;
    handles.accessoric_fig_handle = F;
    S1 = subplot(2,1,1);
    plot(Wavevec(1,:),'r');
    hold on
    plot(Wavevec(2,:),'b');
    plot(Wavevec(3,:),'g');
    plot(Wavevec(4,:),'m');
    legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
    hold off
    y_lim = ylim;
    intlen =  length(Wavevec);
    axis([0 intlen y_lim(1) y_lim(2)]);
    H(1) = get(gcf,'CurrentAxes');
    S2 = subplot(2,1,2);
    handles.subplot2_handle = S2;
    handles.subplot1_handle = S1;
end

% Plot on the accessoric figure - return plot or spike train with structural elements on subplot2
if isequal(index_selected2,1)
    plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected2 + 1;
    plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end
H(2) = get(gcf,'CurrentAxes');
handles.accessoric_axis_handle= H;

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
handles.plot1_handle = plot1_handle;
set(handles.axes1,'ButtonDownFcn',bdf1);

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
if ~ishandle(h)
    h = handles.figure1;
end
guidata(h,handles);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callfftswitch')
set(F,'KeyPressFcn','b_callgreyplot')

% --------------------------------------------------------------------
function varargout = lb_natclust(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected1 = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cl = file_list{index_selected1}(11:30);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
ff = fullfile(pth,fln);
load(ff)
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

d = length(IntraBurstIvVar);
MaxClusterNumber = d;

% Plot on axis1 - variance plot
axes(handles.axes1)
bdf1 = get(handles.axes1,'ButtonDownFcn');
ICC = get(handles.edit1,'String')
Vdisc_new = Vdisc - (i_first - 320000);
b_natclust(Vdisc_new,ICC);

% Create list for listbox2
set(handles.listbox2,'String',[]);

% Plot on axis2 - spike train with structural elements
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
z = zeros(1,length(Time));
z(Vdisc_new) = 1;
z(Vdisc_new+1) = -1;
dec = index_selected + 1;
plot2_handle = plot(Time,z,'b');
axis([Time(1) Time(end) -1.5 1.5])
hold on;
for j = 1:size(Burst{dec},2)
    ind1 = Vdisc_new(Burst{dec}(1,j));
    ind2 = Vdisc_new(Burst{dec}(2,j));
    plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%     rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%     line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%     set(line_handle,'ButtonDownFcn',bdf2)
end

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
handles.plot1_handle = plot1_handle;
set(handles.axes1,'ButtonDownFcn',bdf1);

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
if ~ishandle(h)
    h = handles.figure1;
end
guidata(h,handles);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callfftswitch')
set(F,'KeyPressFcn','b_callgreyplot')



% --------------------------------------------------------------------
% SPIKE TRAIN - RETURN PLOT LISTBOX callback
% --------------------------------------------------------------------
function varargout = listbox2_Callback(h, eventdata, handles, varargin)
switch handles.whichplot
case 'variance'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb2(h, eventdata, handles, varargin)
    end
case 'waveletvector'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb2_wavevec(h, eventdata, handles, varargin)
    end
case 'lombperiodogram'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb2_lombper(h, eventdata, handles, varargin)
    end
case 'fft'
    if ~isempty(find(strcmp(get(handles.figure1,'SelectionType'),{'open','normal'})))
        lb2_fft(h, eventdata, handles, varargin)
    end
end

% --------------------------------------------------------------------
function varargout = lb2(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
if isequal(handles.setinput,['base'])
    nm = '';
else
    nm = ['ica_beta_' handles.setinput];
end
cl = file_list{index_selected}(11:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\',nm];
fln = ['THETA_ICA_' cl];
if isempty(handles.varcell)
    ff = fullfile(pth,fln);
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
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
bdf2 = get(handles.axes2,'ButtonDownFcn');
index_selected = get(handles.listbox2,'Value');
item_list = get(handles.listbox1,'String');
if isequal(index_selected,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
    hans = [];
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc  - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected + 1;
    plot2_handle = plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    sbd = size(Burst{dec},2);
    hans = zeros(1,sbd);
    for j = 1:sbd
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        hans(j) = plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r');
%         rajz = [Time(Vdisc_new(Burst{dec}(1,j))) Time(Vdisc_new(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

% --------------------------------------------------------------------
function varargout = lb2_wavevec(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cl = file_list{index_selected}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

if isempty(handles.varcell)
    ff = fullfile(pth,fln);
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
d = length(IntraBurstIvVar);
MaxClusterNumber = d;

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
index_selected = get(handles.listbox2,'Value');
bdf2 = get(handles.axes2,'ButtonDownFcn');
if isequal(index_selected,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc  - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected + 1;
    plot2_handle = plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown funcions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

% --------------------------------------------------------------------
function varargout = lb2_lombper(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cl = file_list{index_selected}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

if isempty(handles.varcell)
    ff = fullfile(pth,fln);
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
d = length(IntraBurstIvVar);
MaxClusterNumber = d;

% Plot on axis1 - wavelet magnitude average vectors or Lomb periodogram
axes(handles.axes1)
index_selected = get(handles.listbox2,'Value');
bdf1 = get(handles.axes1,'ButtonDownFcn');
if isequal(index_selected,1)
    index_selected_lb1 = get(handles.listbox1,'Value');
    file_list = get(handles.listbox1,'String');
    ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave\',file_list{index_selected_lb1});
    load(ff)
    plot1_handle = plot(Wavevec(1,:),'r');
    hold on
    plot1_handle = plot(Wavevec(2,:),'b');
    plot1_handle = plot(Wavevec(3,:),'g');
    plot1_handle = plot(Wavevec(4,:),'m');
    legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
    hold off
    y_lim = ylim;
    intlen =  length(Wavevec);
    axis([0 intlen y_lim(1) y_lim(2)]);
else
    Vdisc_new = Vdisc - (i_first - 320000);
    dec = index_selected + 1;
    bas = [];
    for bno = 1:size(Burst{dec},2)
        b = Vdisc_new(Burst{dec}(1,bno):Burst{dec}(2,bno));
        bas = [bas b];
    end
%     fbs = Vdisc_new(Burst{dec}(1,:));
%     lenu = length(Time);
%     vdisc_for_lomb = Vdisc_new;
    fname = cl(1:6);
    period_handle = b_lombper_for_ica(fname,i_first,i_second,dec);
    plot1_handle = period_handle;
    handles.plot1_handle = period_handle;
%     [pxx,pyy,jmax,prob,z,effm,period_handle] = b_lomb_period(vdisc_for_lomb,bas,lenu);
end

% Plot on axis2 - return plot or spike train with structural elements
axes(handles.axes2)
index_selected = get(handles.listbox2,'Value');
bdf2 = get(handles.axes2,'ButtonDownFcn');
if isequal(index_selected,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected + 1;
    plot2_handle = plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

ach1 = allchild(handles.axes1);
set(ach1,'ButtonDownFcn',bdf1);
handles.plot1_handle = plot1_handle;
set(handles.axes1,'ButtonDownFcn',bdf1);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callexportfig')

% --------------------------------------------------------------------
function varargout = lb2_fft(h, eventdata, handles, varargin)

% Load ICA data
global DATAPATH
index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
cl = file_list{index_selected}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

if isempty(handles.varcell)
    ff = fullfile(pth,fln);
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
d = length(IntraBurstIvVar);
MaxClusterNumber = d;

% Plot on axis2 - return plot or Lomb periodogram
axes(handles.axes2)
index_selected1 = get(handles.listbox1,'Value');
index_selected2 = get(handles.listbox2,'Value');
bdf2 = get(handles.axes2,'ButtonDownFcn');
if isequal(index_selected2,1)
    plot2_handle = plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    Vdisc_new = Vdisc - (i_first - 320000);
    dec = index_selected2 + 1;
    bas = [];
    for bno = 1:size(Burst{dec},2)
        b = Vdisc_new(Burst{dec}(1,bno):Burst{dec}(2,bno));
        bas = [bas b];
    end
%     fbs = Vdisc_new(Burst{dec}(1,:));
%     lenu = length(Time);
%     vdisc_for_lomb = Vdisc_new;
    fname = cl(1:6);
    period_handle = b_lombper_for_ica(fname,i_first,i_second,dec);
    plot2_handle = period_handle;
%     [pxx,pyy,jmax,prob,z,effm,period_handle] = b_lomb_period(vdisc_for_lomb,bas,lenu);
end

% Plot on the accessoric figure - Wavelet average vectors on subplot1
handles.greyplot1_handle = [];
handles.greyplot2_handle = [];
handles.greyplot3_handle = [];
handles.greyplot4_handle = [];
handles.line_handles = [];
if ~isempty(handles.accessoric_axis_handle)
    H = handles.accessoric_axis_handle;
    F = handles.accessoric_fig_handle;
    S1 = handles.subplot1_handle;
    S2 = handles.subplot2_handle;
    axes(H(2));
else
    F = figure;
    handles.accessoric_fig_handle = F;
    S1 = subplot(2,1,1);
    plot(Wavevec(1,:),'r');
    hold on
    plot(Wavevec(2,:),'b');
    plot(Wavevec(3,:),'g');
    plot(Wavevec(4,:),'m');
    legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
    hold off
    y_lim = ylim;
    intlen =  length(Wavevec);
    axis([0 intlen y_lim(1) y_lim(2)]);
    H(1) = get(gcf,'CurrentAxes');
    S2 = subplot(2,1,2);
    handles.subplot1_handle = S1;
    handles.subplot2_handle = S2;
end

% Plot on the accessoric figure - return plot or spike train with structural elements on subplot2
if isequal(index_selected2,1)
    plot(ReturnPlotXData,ReturnPlotYData,'.');
else
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected2 + 1;
    plot(Time,z,'b');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'r')
%         rajz = [Time(Vdisc(Burst{dec}(1,j))) Time(Vdisc(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','b');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
end
H(2) = get(gcf,'CurrentAxes');
handles.accessoric_axis_handle= H;
guidata(h,handles)

% Set buttondown functions
ach2 = allchild(handles.axes2);
set(ach2,'ButtonDownFcn',bdf2);
handles.plot2_handle = plot2_handle;
set(handles.axes2,'ButtonDownFcn',bdf2);

% Set keypress functions
global HANDLES
HANDLES = handles;
set(handles.figure1,'KeyPressFcn','b_callfftswitch')
set(F,'KeyPressFcn','b_callgreyplot')



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
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% UNDER SUBMENU callback
% --------------------------------------------------------------------
function varargout = under_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'under';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% OVER SUBMENU callback
% --------------------------------------------------------------------
function varargout = over_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'over';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% DECREASE SUBMENU callback
% --------------------------------------------------------------------
function varargout = decrease_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'decrease';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% NODECREASE SUBMENU callback
% --------------------------------------------------------------------
function varargout = nodecrease_submenu_Callback(h, eventdata, handles, varargin)
handles.setinput = 'nodecrease';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin);
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
function varargout = view_menu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
% VARIANCEPLOT SUBMENU callback
% --------------------------------------------------------------------
function varargout = varianceplot_submenu_Callback(h, eventdata, handles, varargin)
handles.whichplot = 'variance';
set(handles.listbox1,'Value',1);
setliststring(h, eventdata, handles, varargin);
lb(h, eventdata, handles, varargin)
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% WAVELETVECTORS SUBMENU callback
% --------------------------------------------------------------------
function varargout = waveletvectors_submenu_Callback(h, eventdata, handles, varargin)
handles.whichplot = 'waveletvector';
set(handles.listbox1,'Value',1)
setliststring_wavevec(h, eventdata, handles, varargin);
lb_wavevec(h, eventdata, handles, varargin)
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% LOMBPERIODOGRAM SUBMENU callback
% --------------------------------------------------------------------
function varargout = lombperiodogram_submenu_Callback(h, eventdata, handles, varargin)
handles.whichplot = 'lombperiodogram';
set(handles.listbox1,'Value',1)
setliststring_wavevec(h, eventdata, handles, varargin);
lb_lombper(h, eventdata, handles, varargin)
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
% FFT SUBMENU callback
% --------------------------------------------------------------------
function varargout = fft_submenu_Callback(h, eventdata, handles, varargin)
handles.whichplot = 'fft';
set(handles.listbox1,'Value',1)
setliststring_wavevec(h, eventdata, handles, varargin);
lb_fft(h, eventdata, handles, varargin)
handles = guidata(h);
guidata(h,handles)



% --------------------------------------------------------------------
function varargout = clustering_submenu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
% NATURAL CLUSTERING SUBMENU callback
% --------------------------------------------------------------------
function varargout = natclust_submenu_Callback(h, eventdata, handles, varargin)
handles.whichplot = 'natclust';
ICC = get(handles.edit1,'String')
ICC = str2num(ICC);
while isempty(ICC) | ICC <= 0 | ICC >=2
    msgbox('Change the cut value of ICC!','Inconsistency Coefficient','warn');
end
lb_natclust(h, eventdata, handles, varargin)
handles = guidata(h);
guidata(h,handles)


% --------------------------------------------------------------------
% ITERATIVE CLUSTERING SUBMENU callback
% --------------------------------------------------------------------
function varargout = itclust_submenu_Callback(h, eventdata, handles, varargin)
handles.whichplot = 'itclust';
lb(h, eventdata, handles, varargin)
handles = guidata(h);
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
        switch handles.whichplot
        case 'variance'
            lb(h, eventdata, handles, varargin)
        case 'waveletvector'
            lb_wavevec(h, eventdata, handles, varargin)
        case 'lombperiodogram'
            lb_lombper(h, eventdata, handles, varargin)
        case 'fft'
            lb_fft(h, eventdata, handles, varargin)
        end
    end
otherwise
    if index_selected - 1 >= 1
        set(handles.listbox1,'Value',index_selected-1)
        switch handles.whichplot
        case 'variance'
            lb(h, eventdata, handles, varargin)
        case 'waveletvector'
            lb_wavevec(h, eventdata, handles, varargin)
        case 'lombperiodogram'
            lb_lombper(h, eventdata, handles, varargin)
        case 'fft'
            lb_fft(h, eventdata, handles, varargin)
        end
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
        switch handles.whichplot
        case 'variance'
            lb2(h, eventdata, handles, varargin)
        case 'waveletvector'
            lb2_wavevec(h, eventdata, handles, varargin)
        case 'lombperiodogram'
            lb2_lombper(h, eventdata, handles, varargin)
        case 'fft'
            lb2_fft(h, eventdata, handles, varargin)
        end
    end
otherwise
    if index_selected - 1 >= 1
        set(handles.listbox2,'Value',index_selected-1)
        switch handles.whichplot
        case 'variance'
            lb2(h, eventdata, handles, varargin)
        case 'waveletvector'
            lb2_wavevec(h, eventdata, handles, varargin)
        case 'lombperiodogram'
            lb2_lombper(h, eventdata, handles, varargin)
        case 'fft'
            lb2_fft(h, eventdata, handles, varargin)
        end
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
ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly',nm);
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



% ---------------------------------------------------------------------
% SETLISTSTRING_WAVEVEC subfunction - set string of listbox1
% ---------------------------------------------------------------------
function varargout = setliststring_wavevec(h, eventdata, handles, varargin)
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave');
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