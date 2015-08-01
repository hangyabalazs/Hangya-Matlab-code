function varargout = b_analyse4(varargin)
% B_ANALYSE4 Application M-file for b_analyse4.fig
%    FIG = B_ANALYSE4 launch b_analyse4 GUI.
%    B_ANALYSE4('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 30-Sep-2009 18:04:34

%ANALYSE   GUI for DSP of EEG and extracellular potential.
%   EEG is displayed on the upper axes, while unit is displayed on the lower
%   one. New registration file can be opened with 'Open' pushbutton. Subsegment
%   can be selected using the slider or editing the text boxes in the lower
%   right corner.
%
%   Help for ANALYSIS submenus:
%       WAVELET
%           EEG - EEG wavelet POWER or PHASE
%           UNIT - UNIT wavelet POWER or PHASE
%           EEG and UNIT - Both EEG and UNIT wavelet POWER or PHASE
%           CROSSWAVELET - crosswavelet power or phase
%           AVERAGE BY SCALE - EEG or UNIT wavelet is split into four submatrices
%               regarding the conventional frequency bands. The submatrices are 
%               averaged along scale and the result is ploted against time. Fourier
%               transformated of these functions can also be computed. You can 
%               switch between frequency bands while watching the Fourier transformated
%               functions using 'q' and 'w' buttons.
%           MODULOGRAM - Fourier transformated is calculated at each frequency of wavelet
%               and the result is displayed on a frequency-against-frequency plot applying
%               a color code where higher peaks of the FFTs are displayed with warmer
%               colors.
%           SETUP - wavelet settings are able to be modified.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       L A U N C H   G U I                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

% Use system color scheme for figure:
set(fig,'Color',get(0,'DefaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
    
    % Store default GUI position in the handles structure for the Resize Function
    old_units = get(handles.figure1,'Units');
    set(handles.figure1,'Units','pixels');
    default_gui_pos = get(handles.figure1,'Position');
    set(handles.figure1,'Units',old_units);
    handles.gui_pos = default_gui_pos;
    
    % Set menus
    ch = get(handles.analysis_menu,'Children');
    set(ch,'Enable','off')
    
    % Create line handle
    set(handles.figure1,'HandleVisibility','on')
    axes(handles.axes_eeg)
    handles.plot1_handle = plot(1,'Color',[1 1 1]);
    axes(handles.axes_unit)
    handles.plot2_handle = plot(1,'Color',[1 1 1]);
    set(handles.figure1,'HandleVisibility','callback')
    
    % Base of logarithm
    handles.la = 8;
    
    % Next Burst button
    handles.burstindex = 0;
    handles.notburstindex = 0;
    handles.burstnum = -1;
    global DATAPATH
    burstdir = [DATAPATH 'LentiFenti\Cluster\mat2\'];
    burstdir0 = [DATAPATH 'LentiFenti\Cluster\mat\'];
    handles.burstdir = burstdir;
    handles.burstdir0 = burstdir0;
    global BURSTDIR
    BURSTDIR = burstdir;
    
    % Set pushbutton 'CData'
    Cstart = imread([matlabroot '\work\Balazs\doubleleft_arrow.bmp']);
    Cprev = imread([matlabroot '\work\Balazs\left_arrow.bmp']);
    Cnext = imread([matlabroot '\work\Balazs\right_arrow.bmp']);
    Clast = imread([matlabroot '\work\Balazs\doubleright_arrow.bmp']);
    set(handles.pushbutton3,'CData',Cstart);
    set(handles.pushbutton5,'CData',Cnext);
    set(handles.pushbutton4,'CData',Cprev);
    set(handles.pushbutton6,'CData',Cstart);
    set(handles.pushbutton8,'CData',Cnext);
    set(handles.pushbutton7,'CData',Cprev);
    
    % Sampling rate
    handles.sr = 20000;
    disp(['Sampling rate: ' num2str(handles.sr)])
    
    % Set default value of whichperiodogram
    global WHICHPERIODOGRAM
    WHICHPERIODOGRAM = 1;
        
    % Load wavelet properties
    ff = fullfile(DATAPATH,'analyse_gui\settings\wavprop');
    try
        load(ff)
        handles.wavelet_properties = wavprop;
        
        default.samprate = 400;
        default.mif = 0.5;
        default.param = 6;
        default.dj = 0.05;
        default.mother = 'morlet';
        default.cfd = 1;
        default.smj1 = 0;
        default.j1 = '';
        default.dst = 1;
        default.sms0 = 0;
        default.s0 = '';
        default.yes = 1;
        default.no = 0;
        % IN CASE YOU MODIFY DEFAULT SETTING, YOU HAVE TO MODIFY THEM IN WAVELET_SETTINGS.M
        % AS WELL!
%         if ~isequal(wavprop,default)
%             msgbox('Wavelet settings differ from default.','Wavelet settings','warn');
%         end
    catch
        msgbox('Cannot load wavelet settings: permission denied or file does not exist. Select Setup from Analysis menu, Wavelet submenu before apply wavelet analysis.',...
            'Wavelet settings','warn');
        newroot = fullfile(DATAPATH,'analyse_gui\');
        settings_dir = fullfile(newroot,'\settings');
        if ~isdir(settings_dir)
            mkdir(newroot,'settings');
        end
    end
    
    % Update 'handles' structure
    guidata(fig,handles);
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
        if ~isempty(lasterr)
            if isempty(findstr(lasterr,'Values must be increasing and non-NaN'))
                disp(lasterr);
                errordlg(lasterr,'Error','modal');
            end
        end
	end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        C A L L B A C K S                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------------------------------------------------------------
% Callback for Resize function
% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
gui_pos = handles.gui_pos;

% Remember the old Units properties
old_units_fig = get(handles.figure1,'Units');
old_units_push = get(handles.pushbutton1,'Units');
old_units_rad = get(handles.radiobutton1,'Units');
old_units_edit = get(handles.edit1,'Units');
old_units_text = get(handles.text1,'Units');
old_units_axes = get(handles.axes_eeg,'Units');
old_units_slider = get(handles.slider1,'Units');
old_units_check = get(handles.checkbox1,'Units');

% Set Units to pixels
set(handles.figure1,'Units','pixels');
pos = get(handles.figure1,'Position');

% Compensate if too narrow
if pos(3) < 150
    if pos(1) == gui_pos(1)
        pos(3) = 150;
        set(handles.figure1,'Position',pos);
    else
        pos(1) = pos(1) + pos(3) - 150;
        pos(3) = 150;
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
for t = 1:8
    eval(['set(handles.pushbutton',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_push{t} = get(handles.pushbutton',int2str(t),',''Position'');']);
    eval(['set(handles.pushbutton',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_push{t}(1) pos(4)-gui_pos(4)+pos_push{t}(2) pos_push{t}(3) pos_push{t}(4)]);']);
    eval(['set(handles.pushbutton',int2str(t),',''Units'',old_units_push);']);
end
for t = 9:10
    eval(['set(handles.pushbutton',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_push{t} = get(handles.pushbutton',int2str(t),',''Position'');']);
    eval(['set(handles.pushbutton',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_push{t}(1) pos_push{t}(2) pos_push{t}(3) pos_push{t}(4)]);']);
    eval(['set(handles.pushbutton',int2str(t),',''Units'',old_units_push);']);
end

% New position of the radiobuttons
pos_rad = cell(1,2);
for t = 1:2
    eval(['set(handles.radiobutton',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_rad{t} = get(handles.radiobutton',int2str(t),',''Position'');']);
    eval(['set(handles.radiobutton',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_rad{t}(1) pos(4)-gui_pos(4)+pos_rad{t}(2) pos_rad{t}(3) pos_rad{t}(4)]);']);
    eval(['set(handles.radiobutton',int2str(t),',''Units'',old_units_rad);']);
end

% New position of the edit texts
pos_edit = cell(1,2);
for t = 1:2
    eval(['set(handles.edit',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_edit{t} = get(handles.edit',int2str(t),',''Position'');']);
    eval(['set(handles.edit',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_edit{t}(1) pos_edit{t}(2) pos_edit{t}(3) pos_edit{t}(4)]);']);
    eval(['set(handles.edit',int2str(t),',''Units'',old_units_edit);']);
end

% New position of the static text boxes above the edit texts and buttons
pos_text = cell(1,5);
for t = 4:5
    eval(['set(handles.text',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_text{t} = get(handles.text',int2str(t),',''Position'');']);
    eval(['set(handles.text',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_text{t}(1) pos(4)-gui_pos(4)+pos_text{t}(2) pos_text{t}(3) pos_text{t}(4)]);']);
    eval(['set(handles.text',int2str(t),',''Units'',old_units_text);']);
end
for t = 3
    eval(['set(handles.text',int2str(t),',''Units'',''pixels'');']);
    eval(['pos_text{t} = get(handles.text',int2str(t),',''Position'');']);
    eval(['set(handles.text',int2str(t),...
            ',''Position'',[pos(3)-gui_pos(3)+pos_text{t}(1) pos_text{t}(2) pos_text{t}(3) pos_text{t}(4)]);']);
    eval(['set(handles.text',int2str(t),',''Units'',old_units_text);']);
end

% New position of the axes
pos_axes = cell(1,2);
set(handles.axes_eeg,'Units','pixels');
pos_axes{1} = get(handles.axes_eeg,'Position');
set(handles.axes_eeg,...
    'Position',[pos_axes{1}(1) pos_axes{1}(2)+((pos(4)-gui_pos(4))/2)...
        pos(3)-gui_pos(3)+pos_axes{1}(3) ((pos(4)-gui_pos(4))/2)+pos_axes{1}(4)]);
set(handles.axes_eeg,'Units',old_units_axes);

set(handles.axes_unit,'Units','pixels');
pos_axes{2} = get(handles.axes_unit,'Position');
set(handles.axes_unit,...
    'Position',[pos_axes{2}(1) pos_axes{2}(2)...
        pos(3)-gui_pos(3)+pos_axes{2}(3) ((pos(4)-gui_pos(4))/2)+pos_axes{2}(4)]);
set(handles.axes_unit,'Units',old_units_axes);

% New position of the static texts above the axes
set(handles.text2,'Units','pixels');
pos_text{2} = get(handles.text2,'Position');
set(handles.text2,...
    'Position',[((pos(3)-gui_pos(3))/2)+pos_text{2}(1) pos(4)-gui_pos(4)+pos_text{2}(2)...
        pos_text{2}(3) pos_text{2}(4)]);
set(handles.text2,'Units',old_units_text);

set(handles.text1,'Units','pixels');
pos_text{1} = get(handles.text1,'Position');
set(handles.text1,...
    'Position',[((pos(3)-gui_pos(3))/2)+pos_text{1}(1) ((pos(4)-gui_pos(4))/2)+pos_text{1}(2)...
        pos_text{1}(3) pos_text{1}(4)]);
set(handles.text1,'Units',old_units_text);

% New position of the slider
set(handles.slider1,'Units','pixels');
pos_slider = get(handles.slider1,'Position');
set(handles.slider1,'Position',...
    [pos_slider(1) pos_slider(2) pos(3)-gui_pos(3)+pos_slider(3) pos_slider(4)]);
set(handles.slider1,'Units',old_units_slider);

% New position of checkbox
set(handles.checkbox1,'Units','pixels');
pos_check = get(handles.checkbox1,'Position');
set(handles.checkbox1,'Position',...
    [pos(3)-gui_pos(3)+pos_check(1) pos_check(2) pos_check(3) pos_check(4)]);
set(handles.checkbox1,'Units',old_units_check);

% Reposition GUI on screen
% movegui(h,'onscreen')

% Reset the Units properties
set(handles.figure1,'Units',old_units_fig);

% Reset gui_pos (GUI position) in the handles structure
handles.gui_pos = pos;
guidata(h,handles);



% --------------------------------------------------------------------
% PUSHBUTTON CALLBACKS
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% Callback for Open button - displays an open dialog
% --------------------------------------------------------------------
function varargout = Open_Callback(h, eventdata, handles, varargin)

% Set directory and open
global DATADIR
global DATAPATH
% cd('X:\In_Vivo\balazs\_analysis\Hajni\EEGMPO\mat\')
% [filename, pathname] = uigetfile({'*.mat', 'All MAT-Files (*.mat)'; ...
% 		'*.*','All Files (*.*)'},'Select Datafile');
cd(DATADIR)
[filename, pathname] = uigetfile({'*.mat', 'All MAT-Files (*.mat)'; ...
    '*.*','All Files (*.*)'},'Select Datafile');

% If "Cancel" is selected then return
% profile on -detail operator
if isequal([filename,pathname],[0,0])
    return
end

% Message
set(handles.figure1,'Name','Loading...');

% Otherwise construct the full filename, load and plot the data
global FILENAME
FILENAME = filename;
fn = fullfile(pathname,filename);
dirstr = [DATAPATH 'analyse_gui\additional_data_for_analyse_gui\'];
nm = [dirstr filename(1:end-4) '_vs'];
dr = dir(dirstr);
cell_dr = cell(1,length(dr)-1);
for m = 3:length(dr),
    cell_dr{m} = dr(m).name(1:end-7);
end
ss = strcmp(cell_dr,filename(1:end-4));
isfs = isempty(find(ss));
if isfs     %first time: import
    [eeg,unit,time] = in_for_analyse(fn,filename,pathname,handles.sr);
    handles.eegs = b_createvs(eeg,handles.la);
    handles.units = b_createvs(unit,handles.la);
    min_unit = min(handles.units{end});
    max_unit = max(handles.units{end});
    min_eeg = min(handles.eegs{end});
    max_eeg = max(handles.eegs{end});
    eeglength = length(eeg);    
else        %if 'vs' already exists: load 'vs'
    load(nm);
    mintafr = handles.sr;
    dt = 1 / mintafr; 
    time = [0:eeglength - 1] * dt; 
    handles.eegs = eegs;
    handles.units = units;
    min_unit = min(units{end});
    max_unit = max(units{end});
    min_eeg = min(eegs{end});
    max_eeg = max(eegs{end});
end

% Plot
axes(handles.axes_eeg);
axis([time(1) time(end) min_eeg max_eeg]);
b_liner(handles.plot1_handle,handles.eegs,handles.gui_pos,handles.la,eeglength,handles.sr);

axes(handles.axes_unit);
axis([time(1) time(end) min_unit max_unit]);
b_liner2(handles.plot2_handle,handles.units,handles.gui_pos,handles.la,eeglength,handles.sr);

drawnow

% Saving 'vs'
if isfs
    eegs = handles.eegs;
    eegs{1} = [];
    units = handles.units;
    units{1} = [];
    eeglength = length(eeg);
    eval(['save(''',nm,''',''eegs'',''units'',''eeglength'')']);
else
    [eeg,unit,time] = in_for_analyse(fn,filename,pathname,handles.sr);
    handles.eegs{1} = eeg;
    handles.units{1} = unit;
end

% Next Burst button
handles.burstindex = 0;
handles.notburstindex = 0;
handles.burstnum = -1;
guidata(h,handles);

% Set slider store variables in 'handles' sturcture
windowwidth = time(end) - time(1);
slidstep = windowwidth / (time(end) - windowwidth + eps - time(1));
set(handles.slider1,'Value',0,'Min',time(1),'Max',time(end)-windowwidth+eps,'SliderStep',...
    [0.01,slidstep],'Enable','off');
handles.fname = filename(1:6);
handles.filename = filename(1:end-4);
handles.pathname = pathname;
handles.time = time;
handles.unit = unit;
handles.eeg = eeg;
handles.len = eeglength;
min_max = [{min_unit} {max_unit} {min_eeg} {max_eeg}];
handles.min_max = min_max;
handles.windowwidth = windowwidth;
guidata(h,handles);
drawnow;

% Set buttondown functions
set_ButtonDownFcns(h, eventdata, handles, varargin)
handles = guidata(h);

% Set edit strings
set_EditStrings(h, eventdata, handles, varargin)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);

% Set menus
ch = get(handles.analysis_menu,'Children');
set(ch,'Enable','on')

% Clear the global variables of previous cell
clear global SCALE WAVETIME WAVPROP_CURRENT HANDLES
clear global WAVE_EEG WAVE_EEG_COORD WAVE_EEG_PROP WAVE_UNIT WAVE_UNIT_COORD WAVE_UNIT_PROP
clear global AVERAGE_EEG AVERAGE_EEG_COORD AVERAGE_UNIT AVERAGE_UNIT_COORD

% Set default value of whichperiodogram
global WHICHPERIODOGRAM
WHICHPERIODOGRAM = 1;

% Clear message
set(handles.figure1,'Name','Raphe Data Analysis');
    
    

% --------------------------------------------------------------------
function [eeg,unit,time] = in_for_analyse(fn,filename,pathname,mintafr)

% Data import
data = load(fn);
if isstruct(data)
    if isfield(data,'EEG')
        eeg = data.EEG.values';
        try
            unit = data.Unit.values';
        catch
            unit = data.unit.values';
        end
    else
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
        if size(data,2) == 1
            data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
        end
        unit = data(:,1)';
        eeg = data(:,2)';
        meret = size(data,1);
    end
else
    if size(data,2) == 1
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
    unit = data(:,1)';
    eeg = data(:,2)';
    meret = size(data,1);
end
clear data
dt = 1 / mintafr;
time = [0:length(unit) - 1] * dt;
xlimit = [min(time),max(time)];

% Filtering unit
disp('Filtering...')
nqf = mintafr / 2;    % Nyquist frequency
flt = fir1(512,100/nqf,'high');
unit = filtfilt(flt,1,unit);



% --------------------------------------------------------------------
% Callback for Discriminate button
% --------------------------------------------------------------------
function varargout = Disc_Callback(h, eventdata, handles, varargin)

% Get threshold
thr = inputdlg('Threshold','Discrimination');
if isempty(thr)
    return
end
thr = thr{1};
if ~isstr(thr)
    return
end
thr = str2num(thr);

% Discriminate
vdisc = disc(handles.unit,thr);
vdisc2 = vdisc;

% Check artefacts
isi = diff(vdisc);
artefact1 = [];
artinx1 = [];
short = find(isi./handles.sr<0.0015);
noiselevel = mean(handles.unit);
if ~isempty(short)
    for k = 1:length(short)
        inx = short(k);
        vdk = vdisc(inx);
        vdkn = vdisc(inx+1);
        un = handles.unit(vdk:vdkn);
        mu = min(un);
        if abs(mu-noiselevel) > abs(mu-(min(handles.unit(vdk),handles.unit(vdkn))))
            vdaa = [vdk vdkn];
            inxaa = [inx inx+1];
            ind = double(handles.unit(vdkn)<handles.unit(vdk)) + 1;
            vda = vdaa(ind);
            inxa = inxaa(ind);
            artefact1 = [artefact1 vda];
            artinx1 = [artinx1 inxa];
        end
    end
end
artefact = unique(artefact1);
artinx = unique(artinx1);
vdisc(artinx1) = [];

% Check big spikes
isi = diff(vdisc);
artefact2 = [];
artinx2 = [];
amp = handles.unit(vdisc);
[samp ix] = sort(amp,'ascend');
next = 1;
while next <= length(amp)
    inx = ix(end-next+1);
    vdk = vdisc(inx);
    setfocus(h,handles,[vdk/handles.sr vdk/handles.sr])
    axes(handles.axes_unit);
    x_lim = get(handles.axes_unit,'XLim');
    y_lim = get(handles.axes_unit,'YLim');
    L = line([x_lim(1)+(x_lim(2)-x_lim(1))/2 x_lim(1)+(x_lim(2)-x_lim(1))/2],...
        [y_lim(2)-0.2 y_lim(2)-0.5],'LineWidth',2,'Color','black');
    bmq = questdlg('Declare as artefact?','Artefact','Yes');
    if isequal(bmq,'Yes')
        artefact2 = [artefact2 vdk];
        artinx2 = [artinx2 inx];
    elseif isequal(bmq,'Cancel')
        delete(L)
        break
    end
    delete(L)
    next = next + 1;
end

% Check small spikes
next = 1;
while next <= length(amp)
    inx = ix(next);
    vdk = vdisc(inx);
    setfocus(h,handles,[vdk/handles.sr vdk/handles.sr])
    axes(handles.axes_unit);
    x_lim = get(handles.axes_unit,'XLim');
    y_lim = get(handles.axes_unit,'YLim');
    L = line([x_lim(1)+(x_lim(2)-x_lim(1))/2 x_lim(1)+(x_lim(2)-x_lim(1))/2],...
        [y_lim(2)-0.2 y_lim(2)-0.5],'LineWidth',2,'Color','black');
    smq = questdlg('Declare as noise?','Noise','Yes');
    if isequal(smq,'Yes')
        artefact2 = [artefact2 vdk];
        artinx2 = [artinx2 inx];
    elseif isequal(smq,'Cancel')
        delete(L)
        break
    end
    delete(L)
    next = next + 1;
end
artefact2 = unique(artefact2);
artinx2 = unique(artinx2);
vdisc(artinx2) = [];

% Check small not-spikes
vdisclow = disc(handles.unit,thr*0.8);
vdiscdiff = setdiff(vdisclow,vdisc2);
adds = [];
if ~isempty(vdiscdiff)
    amp = handles.unit(vdiscdiff);
    [samp ix] = sort(amp,'descend');
    next = 1;
    while next <= length(amp)
        inx = ix(next);
        vdk = vdiscdiff(inx);
        setfocus(h,handles,[vdk/handles.sr vdk/handles.sr])
        axes(handles.axes_unit);
        x_lim = get(handles.axes_unit,'XLim');
        y_lim = get(handles.axes_unit,'YLim');
        L = line([x_lim(1)+(x_lim(2)-x_lim(1))/2 x_lim(1)+(x_lim(2)-x_lim(1))/2],...
            [y_lim(2)-0.2 y_lim(2)-0.5],'LineWidth',2,'Color','black');
        smq = questdlg('Declare as spike?','Noise','Yes');
        if isequal(smq,'Yes')
            adds = [adds vdk];
        elseif isequal(smq,'Cancel')
            delete(L)
            break
        end
        delete(L)
        next = next + 1;
    end
end
vdisc = [vdisc adds];
vdisc = sort(vdisc,'ascend');

% Check artefacts again
isi = diff(vdisc);
artefact3 = [];
artinx3 = [];
short = find(isi./handles.sr<0.0015);
noiselevel = mean(handles.unit);
if ~isempty(short)
    for k = 1:length(short)
        inx = short(k);
        vdk = vdisc(inx);
        vdkn = vdisc(inx+1);
        un = handles.unit(vdk:vdkn);
        mu = min(un);
        if abs(mu-noiselevel) > abs(mu-(min(handles.unit(vdk),handles.unit(vdkn))))
            vdaa = [vdk vdkn];
            inxaa = [inx inx+1];
            ind = double(handles.unit(vdkn)<handles.unit(vdk)) + 1;
            vda = vdaa(ind);
            inxa = inxaa(ind);
            artefact3 = [artefact3 vda];
            artinx3 = [artinx3 inxa];
        end
    end
end
artefact3 = unique(artefact3);
artinx3 = unique(artinx3);
vdisc(artinx3) = [];
artefact = sort(unique([artefact1 artefact2 artefact2]));

% Save
nm = [handles.filename '_d.mat'];
mm = pwd;
cd(handles.pathname)
str = ['save ' nm ' vdisc thr artefact adds'];
eval(str)
cd(mm)



% --------------------------------------------------------------------
% Callback for Prev/Next Burst button
% --------------------------------------------------------------------
function varargout = bnext_Callback(h, eventdata, handles, varargin)

% Set burst index
if ~isequal(handles.burstnum,-1)
    handles.burstindex = min(handles.burstindex+1,handles.burstnum);
else
    handles.burstindex = handles.burstindex + 1;
end
guidata(h,handles)

% Next/prev burst
bprevnext(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = bprev_Callback(h, eventdata, handles, varargin)

% Set burst index
handles.burstindex = max(handles.burstindex-1,1);
guidata(h,handles)

% Next/prev burst
bprevnext(h, eventdata, handles, varargin)



% -------------------------------------------------------------------------
function bprevnext(h, eventdata, handles, varargin)

% Load burst data
global FILENAME
fn = FILENAME;
ff = fullfile(handles.burstdir,[fn(1:end-4) '_CLUST2.mat']);
if exist(ff,'file')
    load(ff)
else
    return
end
vburst = vdisc(Burst);
if isequal(size(vburst),[1,2])
    vburst = vburst';
end
burst = Burst;

% Get current burst
v = get(handles.radiobutton1,'Value');
burstnum = size(vburst,2);
if v
    cvburst = vburst(:,handles.burstindex);
else
    intraburstiv = [];
    firstiv = [];
    for j = 1:burstnum    %computing intraburstiv and firstiv
        b = vdisc(burst(1,j):burst(2,j));
        db = diff(b);
        intraburstiv = [intraburstiv db];
        firstiv = [firstiv db(1)];
    end
    firstiv2 = firstiv + linspace(0.001,0.002,length(firstiv)); % avoid intervals of the same length
    sfi = sort(firstiv2,2,'descend');
    cfi = sfi(handles.burstindex);
    inx = find(firstiv2==cfi);
    cvburst = vburst(:,inx);
end
handles.burstnum = burstnum;
guidata(h,handles)

% Set focus
setfocus(h,handles,[cvburst(1)/handles.sr cvburst(2)/handles.sr])



% --------------------------------------------------------------------
% Callback for Prev/Next Notburst button
% --------------------------------------------------------------------
function varargout = nbnext_Callback(h, eventdata, handles, varargin)

% Set burst index
handles.notburstindex = handles.notburstindex + 1;
guidata(h,handles)

% Next/prev burst
nbprevnext(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = nbprev_Callback(h, eventdata, handles, varargin)

% Set burst index
handles.notburstindex = max(handles.notburstindex-1,1);
guidata(h,handles)

% Next/prev burst
nbprevnext(h, eventdata, handles, varargin)



% -------------------------------------------------------------------------
function nbprevnext(h, eventdata, handles, varargin)

% Load burst data
global FILENAME
fn = FILENAME;
ff = fullfile(handles.burstdir,[fn(1:end-4) '_CLUST2.mat']);
if exist(ff,'file')
    load(ff)
    vburst = vdisc(Burst);
    if isequal(size(vburst),[1,2])
        vburst = vburst';
    end
    burst = Burst;
    isi = diff(vdisc);

    burstnum = size(vburst,2);
    notburstiv = [];
    notburst = [];
    if ~isequal(burst(1,1),1)
        nb = vdisc(1:burst(1,1));
        dnb = diff(nb);
        notburstiv = [notburstiv dnb];
        notburst = [notburst nb(1:end-1)];
    end
    for j = 1:burstnum-1    %computing intraburstiv and firstiv
        nb = vdisc(burst(2,j):burst(1,j+1));
        dnb = diff(nb);
        notburstiv = [notburstiv dnb];
        notburst = [notburst nb(1:end-1)];
    end
    if ~isequal(burst(2,end),length(vdisc))
        nb = vdisc(burst(2,end):end);
        dnb = diff(nb);
        notburstiv = [notburstiv dnb];
        notburst = [notburst nb(1:end-1)];
    end
else
    ff = fullfile(handles.burstdir0,[fn(1:end-4) '_CLUST.mat']);
    load(ff)
    burstnum = 0;
    notburstiv = diff(vdisc);
    notburst = vdisc(1:end-1);
end
notburstiv2 = notburstiv + linspace(0.001,0.002,length(notburstiv)); % avoid intervals of the same length
snbi = sort(notburstiv2,2,'ascend');
cnbi = snbi(handles.notburstindex);
inx = find(notburstiv2==cnbi);
cvnotburst = [notburst(inx) notburst(inx+1)];
handles.burstnum = burstnum;
guidata(h,handles)

% Set focus
setfocus(h,handles,[cvnotburst(1)/handles.sr cvnotburst(2)/handles.sr])



% --------------------------------------------------------------------
% Callback for reset buttons
% --------------------------------------------------------------------
function bstart_Callback(h, eventdata, handles, varargin)

% Reset burstindex
handles.burstindex = 1;
guidata(h,handles)

% Next/prev burst
bprevnext(h, eventdata, handles, varargin)



% -------------------------------------------------------------------------
function nbstart_Callback(h, eventdata, handles, varargin)

% Reset notburstindex
handles.notburstindex = 1;
guidata(h,handles)

% Next/prev notburst
nbprevnext(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
% Callback for Close Figs button - close all
% --------------------------------------------------------------------
function varargout = Close_Figs_Callback(h, eventdata, handles, varargin)
default = get(handles.figure1,'HandleVisibility');
set(handles.figure1,'HandleVisibility','off')
close all
set(handles.figure1,'HandleVisibility',default)



% --------------------------------------------------------------------
% Callback for Close button and CloseRequestFcn - deletes GUI
% --------------------------------------------------------------------
function varargout = Close_Callback(h, eventdata, handles, varargin)
clear global SCALE WAVETIME WHICHPERIODOGRAM WAVPROP_CURRENT HANDLES
clear global WAVE_EEG WAVE_EEG_COORD WAVE_EEG_PROP WAVE_UNIT WAVE_UNIT_COORD WAVE_UNIT_PROP
clear global AVERAGE_EEG AVERAGE_EEG_COORD AVERAGE_UNIT AVERAGE_UNIT_COORD
clear global BURSTDIR FILENAME
delete(handles.figure1)



% --------------------------------------------------------------------
% RADIOBUTTON CALLBACKS
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% Callback for radiobuttons associated to Next Burst button
% --------------------------------------------------------------------
function radiobutton1_Callback(hObject, eventdata, handles)

set(handles.radiobutton1,'Value',1);
set(handles.radiobutton2,'Value',0);



% -------------------------------------------------------------------------
function radiobutton2_Callback(hObject, eventdata, handles)

set(handles.radiobutton2,'Value',1);
set(handles.radiobutton1,'Value',0);



% --------------------------------------------------------------------
% EDIT STRING CALLBACKS
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% Callback for Edit string - set interval's lower limit
% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)

% Get input and limits
x1 = get(handles.edit1,'String');
x2 = get(handles.edit2,'String');
x_lim = get(handles.axes_eeg,'XLim');
x_limr = round(x_lim*10) / 10;
x1 = str2num(x1);
x2 = str2num(x2);
if isempty(x1)
    errordlg('Input must be numeric.','Error','modal');
    set(handles.edit1,'String',num2str(x_limr(1)));
    return
end
if x1 >= x2
    errordlg('Input values must be increasing.','Error','modal');
    set(handles.edit1,'String',num2str(x_limr(1)));
    return
end
y_eeg = get(handles.axes_eeg,'YLim');
y_unit = get(handles.axes_unit,'YLim');

% If out of range, set string back; else set axis and update plot
if x1 < 0 | x1 > handles.len
    set(handles.edit1,'String',num2str(x_limr(1)));
else
    axes(handles.axes_eeg)
    axis([x1 x_lim(2) y_eeg(1) y_eeg(2)]);
    b_liner(handles.plot1_handle,handles.eegs,handles.gui_pos,handles.la,handles.len,handles.sr);
    
    axes(handles.axes_unit)
    axis([x1 x_lim(2) y_unit(1) y_unit(2)]);
    b_liner2(handles.plot2_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
end
drawnow

% Set slider
set_Slider(h, eventdata, handles, varargin)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);



% --------------------------------------------------------------------
% Callback for Edit string - set interval's upper limit
% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)

% Get input and limits
x1 = get(handles.edit1,'String');
x2 = get(handles.edit2,'String');
x_lim = get(handles.axes_eeg,'XLim');
x_limr = round(x_lim*10) / 10;
x2 = str2num(x2);
x1 = str2num(x1);
if isempty(x2)
    errordlg('Input must be numeric.','Error','modal');
    set(handles.edit2,'String',num2str(x_limr(2)));
    return
end
if x1 >= x2
    errordlg('Input values must be increasing.','Error','modal');
    set(handles.edit2,'String',num2str(x_limr(2)));
    return
end
y_eeg = get(handles.axes_eeg,'YLim');
y_unit = get(handles.axes_unit,'YLim');

% If out of range, set string back; else set axis and update plot
if x2 < 0 | x2 > handles.len
    set(handles.edit2,'String',num2str(x_limr(2)));
else
    axes(handles.axes_eeg)
    axis([x_lim(1) x2 y_eeg(1) y_eeg(2)]);
    b_liner(handles.plot1_handle,handles.eegs,handles.gui_pos,handles.la,handles.len,handles.sr);
    
    axes(handles.axes_unit)
    axis([x_lim(1) x2 y_unit(1) y_unit(2)]);
    b_liner2(handles.plot2_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);
end
drawnow

% Set slider
set_Slider(h, eventdata, handles, varargin)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);



% --------------------------------------------------------------------
% SLIDER CALLBACK
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% Callback for Slider
% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)

% Get variables from 'handles' structure
time = handles.time;
eeg = handles.eeg;
unit = handles.unit;
min_max = handles.min_max;
min_unit = min_max{1};
max_unit = min_max{2};
min_eeg = min_max{3};
max_eeg = min_max{4};    
windowwidth = handles.windowwidth;
slider_time = get(h,'Value');

% Make sure that the selected area is within bounds of the data.
% Extract the local x- and y-data for the axes
if slider_time == 0
    slider_time = eps;
end
y_eeg = get(handles.axes_eeg,'YLim');
y_unit = get(handles.axes_unit,'YLim');

% Plot
axes(handles.axes_eeg);
axis([slider_time slider_time+windowwidth y_eeg(1) y_eeg(2)]);
b_liner(handles.plot1_handle,handles.eegs,handles.gui_pos,handles.la,handles.len,handles.sr);

axes(handles.axes_unit);
axis([slider_time slider_time+windowwidth y_unit(1) y_unit(2)]);
b_liner2(handles.plot2_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);

drawnow;

% Set edit strings
set_EditStrings(h, eventdata, handles, varargin)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);



% --------------------------------------------------------------------
% KEYPRESS AND BUTTON DOWN FUNCTIONS
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% Callback for Axes ButtonDownFunction - ZOOMING
% --------------------------------------------------------------------
function varargout = axes_eeg_ButtonDownFcn(h, eventdata, handles, varargin)

% Get variables from 'handles' structure
time = handles.time;
min_max = handles.min_max;
min_unit = min_max{1};
max_unit = min_max{2};
min_eeg = min_max{3};
max_eeg = min_max{4}; 

% Set axis
seltyp = get(handles.figure1,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(handles.axes_eeg,'CurrentPoint'); % button down detected
    units = get(handles.figure1,'units');
    set(handles.figure1,'units','pixels')
    rbbox([get(handles.figure1,'currentpoint') 0 0],get(handles.figure1,'currentpoint'),handles.figure1);                   % return figure units
    set(handles.figure1,'units',units)
    point2 = get(handles.axes_eeg,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    if isequal(point1,point2)
        xx = get(handles.axes_eeg,'XLim');
        yy = get(handles.axes_eeg,'YLim');
        xx2 = (abs(xx(2) - xx(1))) / 4;
        yy2 = (abs(yy(2) - yy(1))) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        if time(1) < xx3(1) & time(end) > xx3(2)
            xlim_unit = xx3;
            xlim_eeg = xx3;
        elseif time(1) > xx3(1)
            xx_new(1) = time(1);
            xx_new(2) = time(1) + (2 * xx2);
            xlim_unit = xx_new;
            xlim_eeg = xx_new;
        elseif time(end) < xx3(2)
            xx_new(1) = time(end) - (2 * xx2);
            xx_new(2) = time(end);
            xlim_unit = xx_new;
            xlim_eeg = xx_new;
        end
        if min_eeg < yy3(1) & max_eeg > yy3(2)
            ylim_unit = yy3;
            ylim_eeg = yy3;
        elseif min_eeg > yy3(1)
            yy_new(1) = min_eeg;
            yy_new(2) = min_eeg + (2 * yy2);
            ylim_unit = yy_new;
            ylim_eeg = yy_new;
        elseif max_eeg < yy3(2)
            yy_new(1) = max_eeg - (2 * yy2);
            yy_new(2) = max_eeg;
            ylim_unit = yy_new;
            ylim_eeg = yy_new;
        end
    else
        xlim_unit = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        xlim_eeg = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        ylim_unit = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
        ylim_eeg = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
    end

case 'open'   % set default
    axes(handles.axes_eeg);
    xlim_unit = [time(1) time(end)];
    xlim_eeg = [time(1) time(end)];
    ylim_unit = [min_unit max_unit];
    ylim_eeg = [min_eeg max_eeg];
    
otherwise   % zoom out
    point = get(handles.axes_eeg,'CurrentPoint'); % button down detected
    point = point(1,1:2);
    xx = get(handles.axes_eeg,'XLim');
    yy = get(handles.axes_eeg,'YLim');
    xx2 = abs(xx(2) - xx(1));
    yy2 = abs(yy(2) - yy(1));
    if xx2 > (time(end) - time(1)) / 2,
        xx2 = (time(end) - time(1)) / 2;
    end
    if yy2 > (abs(max_eeg - min_eeg)) / 2,
        yy2 = (abs(max_eeg - min_eeg)) / 2;
    end
    xx3(1) = point(1) - xx2;
    xx3(2) = point(1) + xx2;
    yy3(1) = point(2) - yy2;
    yy3(2) = point(2) + yy2;
    if time(1) < xx3(1) & time(end) > xx3(2)
        xlim_unit = xx3;
        xlim_eeg = xx3;
    elseif time(1) > xx3(1)
        xx_new(1) = time(1);
        xx_new(2) = time(1) + (2 * xx2);
        xlim_unit = xx_new;
        xlim_eeg = xx_new;
    elseif time(end) < xx3(2)
        xx_new(1) = time(end) - (2 * xx2);
        xx_new(2) = time(end);
        xlim_unit = xx_new;
        xlim_eeg = xx_new;
    end
    if min_eeg < yy3(1) & max_eeg > yy3(2)
        ylim_unit = yy3;
        ylim_eeg = yy3;
    elseif min_eeg > yy3(1)
        yy_new(1) = min_eeg;
        yy_new(2) = min_eeg + (2 * yy2);
        ylim_unit = yy_new;
        ylim_eeg = yy_new;
    elseif max_eeg < yy3(2)
        yy_new(1) = max_eeg - (2 * yy2);
        yy_new(2) = max_eeg;
        ylim_unit = yy_new;
        ylim_eeg = yy_new;
    end
end
set(handles.axes_eeg,'XLim',xlim_eeg);
set(handles.axes_eeg,'YLim',ylim_eeg);
set(handles.axes_unit,'XLim',xlim_unit);
if isequal(get(handles.checkbox1,'Value'),0)
    set(handles.axes_unit,'YLim',ylim_unit);
end
axes(handles.axes_eeg);
b_liner(handles.plot1_handle, handles.eegs,handles.gui_pos,handles.la,handles.len,handles.sr);
axes(handles.axes_unit);
b_liner2(handles.plot2_handle, handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);

exl = get(handles.axes_eeg,'XLim');
uxl = get(handles.axes_unit,'XLim');
set(handles.axes_eeg,'XLim',min(max(exl,time(1)),time(end)));
set(handles.axes_unit,'XLim',min(max(uxl,time(1)),time(end)));
drawnow;

% Set slider
set_Slider(h, eventdata, handles, varargin)
handles = guidata(h);

% Set edit strings
set_EditStrings(h, eventdata, handles, varargin)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);



% --------------------------------------------------------------------
% Callback for Axes ButtonDownFunction - ZOOMING
% --------------------------------------------------------------------
function varargout = axes_unit_ButtonDownFcn(h, eventdata, handles, varargin)

% Get variables from 'handles' structure
time = handles.time;
min_max = handles.min_max;
min_unit = min_max{1};
max_unit = min_max{2};
min_eeg = min_max{3};
max_eeg = min_max{4}; 

% Set axis
seltyp = get(handles.figure1,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(handles.axes_unit,'CurrentPoint'); % button down detected
    units = get(handles.figure1,'units');
    set(handles.figure1,'units','pixels')
    rbbox([get(handles.figure1,'currentpoint') 0 0],get(handles.figure1,'currentpoint'),handles.figure1);                   % return figure units
    set(handles.figure1,'units',units)
    point2 = get(handles.axes_unit,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    if isequal(point1,point2)
        xx = get(handles.axes_unit,'XLim');
        yy = get(handles.axes_unit,'YLim');
        xx2 = (abs(xx(2) - xx(1))) / 4;
        yy2 = (abs(yy(2) - yy(1))) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        if time(1) < xx3(1) & time(end) > xx3(2)
            xlim_unit = xx3;
            xlim_eeg = xx3;
        elseif time(1) > xx3(1)
            xx_new(1) = time(1);
            xx_new(2) = time(1) + (2 * xx2);
            xlim_unit = xx_new;
            xlim_eeg = xx_new;
        elseif time(end) < xx3(2)
            xx_new(1) = time(end) - (2 * xx2);
            xx_new(2) = time(end);
            xlim_unit = xx_new;
            xlim_eeg = xx_new;
        end
        if min_eeg < yy3(1) & max_eeg > yy3(2)
            ylim_unit = yy3;
            ylim_eeg = yy3;
        elseif min_eeg > yy3(1)
            yy_new(1) = min_eeg;
            yy_new(2) = min_eeg + (2 * yy2);
            ylim_unit = yy_new;
            ylim_eeg = yy_new;
        elseif max_eeg < yy3(2)
            yy_new(1) = max_eeg - (2 * yy2);
            yy_new(2) = max_eeg;
            ylim_unit = yy_new;
            ylim_eeg = yy_new;
        end
    else
        xlim_unit = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        xlim_eeg = [min([point1(1) point2(1)]) max([point1(1) point2(1)])];
        ylim_unit = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
        ylim_eeg = [min([point1(2) point2(2)]) max([point1(2) point2(2)])];
    end

case 'open'   % set default
    axes(handles.axes_eeg);
    xlim_unit = [time(1) time(end)];
    xlim_eeg = [time(1) time(end)];
    ylim_unit = [min_unit max_unit];
    ylim_eeg = [min_eeg max_eeg];
    
otherwise   % zoom out
    point = get(handles.axes_unit,'CurrentPoint'); % button down detected
    point = point(1,1:2);
    xx = get(handles.axes_unit,'XLim');
    yy = get(handles.axes_unit,'YLim');
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
        xlim_unit = xx3;
        xlim_eeg = xx3;
    elseif time(1) > xx3(1)
        xx_new(1) = time(1);
        xx_new(2) = time(1) + (2 * xx2);
        xlim_unit = xx_new;
        xlim_eeg = xx_new;
    elseif time(end) < xx3(2)
        xx_new(1) = time(end) - (2 * xx2);
        xx_new(2) = time(end);
        xlim_unit = xx_new;
        xlim_eeg = xx_new;
    end
    if min_unit < yy3(1) & max_unit > yy3(2)
        ylim_unit = yy3;
        ylim_eeg = yy3;
    elseif min_unit > yy3(1)
        yy_new(1) = min_unit;
        yy_new(2) = min_unit + (2 * yy2);
        ylim_unit = yy_new;
        ylim_eeg = yy_new;
    elseif max_unit < yy3(2)
        yy_new(1) = max_unit - (2 * yy2);
        yy_new(2) = max_unit;
        ylim_unit = yy_new;
        ylim_eeg = yy_new;
    end
end
set(handles.axes_eeg,'XLim',xlim_eeg);
if isequal(get(handles.checkbox1,'Value'),0)
    set(handles.axes_eeg,'YLim',ylim_eeg);
end
set(handles.axes_unit,'XLim',xlim_unit);
set(handles.axes_unit,'YLim',ylim_unit);
axes(handles.axes_eeg);
b_liner(handles.plot1_handle, handles.eegs,handles.gui_pos,handles.la,handles.len,handles.sr);
axes(handles.axes_unit);
b_liner2(handles.plot2_handle, handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);

exl = get(handles.axes_eeg,'XLim');
uxl = get(handles.axes_unit,'XLim');
set(handles.axes_eeg,'XLim',min(max(exl,time(1)),time(end)));
set(handles.axes_unit,'XLim',min(max(uxl,time(1)),time(end)));
drawnow;

% Set slider
set_Slider(h, eventdata, handles, varargin)
handles = guidata(h);

% Set edit strings
set_EditStrings(h, eventdata, handles, varargin)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);



% --------------------------------------------------------------------
% Average fft KeypressDownFunction - switch between frequency bands
% --------------------------------------------------------------------
function b_callfftswitch(st)
% Response to keypress:
%     q - switches up one wavelet average magnitude fft plot
%     w - switches down one wavelet average magnitude fft plot
%     z - switches 'zoom xon' on

inp = get(gcf,'CurrentCharacter');
switch inp
case 'q'
    b_fftswitch('u',st)
case 'w'
    b_fftswitch('d',st)
case 'z'
    zoom xon
end

% --------------------------------------------------------------------
function b_fftswitch(s,st)

% Switch
global WHICHPERIODOGRAM
whichperiod = WHICHPERIODOGRAM;
switch s
case 'u'
    whichperiod = whichperiod + 1;
    whichperiod = mod(whichperiod-1,4) + 1;
case 'd'
    whichperiod = whichperiod - 1;
    whichperiod = mod(whichperiod-1,4) + 1;
end
WHICHPERIODOGRAM = whichperiod;
global WHICHPERIODOGRAM

% Get wavelet average
if isequal(st,'eeg')
    global AVERAGE_EEG
    average = AVERAGE_EEG;
else
    global AVERAGE_UNIT
    average = AVERAGE_UNIT;
end
lenwv = size(average,2);

% Get 'handles' structure
global HANDLES
handles = HANDLES;

% Sapmling rate
wavprop = handles.wavelet_properties;
samprate = wavprop.samprate;

% Plot fft
switch whichperiod
case 1
    [p_1 f_1] = b_fft2(average(4,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it1 - 3 Hz}')
case 2
    [p_2 f_2] = b_fft2(average(3,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it3 - 6 Hz}')
case 3
    [p_3 f_3] = b_fft2(average(2,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it6 - 20 Hz}')
case 4
    [p_4 f_4] = b_fft2(average(1,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it20 - 50 Hz}')
end
tt = getappdata(gcf,'title');
title(tt)



% --------------------------------------------------------------------
% MENU CALLBACKS - WAVELET
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% EEG submenu callbacks
% --------------------------------------------------------------------
function varargout = eegpower_submenu_Callback(h, eventdata, handles, varargin)
wavelet_for_analyse(h,handles,'eeg','power');

% --------------------------------------------------------------------
function varargout = eegphase_submenu_Callback(h, eventdata, handles, varargin)
wavelet_for_analyse(h,handles,'eeg','phase');



% --------------------------------------------------------------------
% Unit submenu callbacks
% --------------------------------------------------------------------
function varargout = unitpower_submenu_Callback(h, eventdata, handles, varargin)
wavelet_for_analyse(h,handles,'unit','power');

% --------------------------------------------------------------------
function varargout = unitphase_submenu_Callback(h, eventdata, handles, varargin)
wavelet_for_analyse(h,handles,'unit','phase');



% --------------------------------------------------------------------
% EEG and Unit ('both') submenu callbacks
% --------------------------------------------------------------------
function varargout = bothpower_submenu_Callback(h, eventdata, handles, varargin)
wavelet_for_analyse(h,handles,'eeg','power');
wavelet_for_analyse(h,handles,'unit','power');

% --------------------------------------------------------------------
function varargout = bothphase_submenu_Callback(h, eventdata, handles, varargin)
wavelet_for_analyse(h,handles,'eeg','phase');
wavelet_for_analyse(h,handles,'unit','phase');



% --------------------------------------------------------------------
% WAVELET function
% --------------------------------------------------------------------
function wavelet_for_analyse(h,handles,st,porp)

% Get wavelet settings
global DATAPATH
try
    ff = fullfile(DATAPATH,'analyse_gui\settings\wavprop');
    load(ff)
    handles.wavelet_properties = wavprop;
    guidata(h,handles)
catch
    msgbox('Cannot load wavelet settings: permission denied or file does not exist. Select Setup from Analysis menu, Wavelet submenu before apply wavelet analysis.',...
        'Wavelet settings','warn');
    return
end

% Create variables for wavelet
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave = wavecheck_in_wavelet(st,datinx1,datinx2);

% Create scale vector
global SCALE        % if 'f' is empty or an inadecvate 'f' is imported,
f = SCALE;          % 'f' is overwritten at the end of wavelet calculation
global WAVETIME     % the same comment for 'wavetime'
wavetime = WAVETIME;

% Wavelet calculation
if isempty(wave)
    str = ['sst = handles.' st '(datinx1:datinx2);'];
    eval(str);
    newstep = handles.sr / wavprop.samprate;
    lensst = length(sst);
    
    if strcmp(st,'unit')        % unit sinc convolution
%         sst = b_sincconv3(sst,newstep,handles.sr);
        sst = sst(1:newstep:lensst);    % downsample eeg
    else
        sst = sst(1:newstep:lensst);    % downsample eeg
    end
    
    sst = (sst - mean(sst)) / std(sst); % standardization
    
    dt = 1 / wavprop.samprate;          % wavelet transformation
    wavetime = [0:length(sst)-1] * dt;
    n = length(sst);
    if wavprop.yes
        pad = 1;
    else
        pad = 0;
    end
    dj = wavprop.dj;
    if wavprop.dst
        s0 = 2 * dt;
    else
        s0 = wavprop.s0;
    end
    if wavprop.cfd
        j1 = ((1 / dj) * log2(n/2)) * 2;
    else
        j1 = wavprop.j1;
    end
    j1 = ceil(j1);
    j = (0:j1);
    s = s0 .* 2 .^ (j * dj);
    omega0 = 6;
    c = 4 * pi / (omega0 + sqrt(2+omega0^2));
    fperiod = c .* s;
    f = 1 ./ fperiod;
    lag1 = 0.72;  
    mother = wavprop.mother;
    param = wavprop.param;
    mif = wavprop.mif;          %minimal intresting frequency
    mis = find(f>mif);
    mis = mis(end);     %maximal intristing scale
    [wave,period,scale,coi] = b_wavelet_new3(sst,dt,pad,dj,s0,j1,mother,param,mis);
end

% Plot
if porp == 'power'
    im = abs(wave).^2;
else
    im = angle(wave);
end
figure
imagesc(im)
ff = round(f*100) / 100;
time = round(wavetime*100) / 100;
b_rescaleaxis('Y',ff)
b_rescaleaxis('X',time)
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' ' upper(st) ' Wavelet'];
title(tt)

drawnow

% Export 'f' and 'time' for ZOOM_FOR_WAVELET and others
global SCALE
SCALE = f;
global WAVETIME
WAVETIME = wavetime;
setappdata(gca,'scalex',wavetime)
setappdata(gca,'scaley',ff)
b_zoomset_for_wavelet

% Export wavelet, wavelet coordinates and wavelet properties
nm = ['wave_' st];
nm = upper(nm);
nmcoord = [nm '_coord'];
nmcoord = upper(nmcoord);
nmprop = [nm '_prop'];
nmprop = upper(nmprop);
eval(['global ' nm]);
eval([nm ' =  wave;']);
eval(['global ' nmcoord]);
eval([nmcoord ' =  [datinx1 datinx2];']);
eval(['global ' nmprop]);
eval([nmprop ' =  wavprop;']);

% --------------------------------------------------------------------
function wave = wavecheck_in_wavelet(st,datinx1,datinx2)
nm = ['wave_' st];
nm = upper(nm);
eval(['global ' nm]);
eval(['wave = ' nm ';']);
if ~isempty(wave)
    nmcoord = [nm '_coord'];
    nmcoord = upper(nmcoord);
    eval(['global ' nmcoord]);
    eval(['coord = ' nmcoord ';']);
    if ~isequal([datinx1 datinx2],coord)
        wave = [];
    else
        nmprop = [nm '_prop'];
        nmprop = upper(nmprop);
        eval(['global ' nmprop]);
        eval(['prop = ' nmprop ';']);
        global WAVPROP_CURRENT
        newprop = WAVPROP_CURRENT;
        if ~isequal(prop,newprop) & ~isempty(newprop)
            wave = [];
        end
    end
end



% --------------------------------------------------------------------
% Crosswavelet submenu callbacks
% --------------------------------------------------------------------
function varargout = crosspower_submenu_Callback(h, eventdata, handles, varargin)
crosswavelet_for_analyse(h,eventdata,handles,varargin,'power')

% --------------------------------------------------------------------
function varargout = crossphase_submenu_Callback(h, eventdata, handles, varargin)
crosswavelet_for_analyse(h,eventdata,handles,varargin,'phase')



% --------------------------------------------------------------------
% CROSSWAVELET function
% --------------------------------------------------------------------
function crosswavelet_for_analyse(h,eventdata,handles,varargin,porp)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_eeg = get_eeg_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);
wave_unit = get_unit_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Calculate crosswavelet
wave_cross = wave_eeg .* conj(wave_unit);

% Plot
if porp == 'power'
    im = abs(wave_cross).^2;
else
    im = angle(wave_cross);
end
H = figure;
imagesc(im)
mx = max(im(:));
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
text(xcoord,ycoord,num2str(mx),'Color',[1 1 1]);
global SCALE    % import 'f' and 'time'
f = SCALE;
global WAVETIME
wavetime = WAVETIME;
ff = round(f*100) / 100;
time = round(wavetime*100) / 100;
b_rescaleaxis('Y',ff)
b_rescaleaxis('X',time)
tt = handles.filename;
tt(tt=='_') = ' ';
tt = [tt ' Crosswavelet ' porp];
title(tt)
setappdata(gca,'scalex',wavetime)
setappdata(gca,'scaley',f)
b_zoomset_for_wavelet

drawnow

% dn = 'X:\In_Vivo\balazs\_analysis\Hajni\EEGMPO\crosswavelet\';
% fn = [dn tt '.fig'];
% saveas(H,fn)



% --------------------------------------------------------------------
% Wavelet average submenu callbacks
% --------------------------------------------------------------------
function varargout = eeg_average_submenu_Callback(h, eventdata, handles, varargin)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_eeg = get_eeg_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Calculate averages and plot
wavelet_average(handles,wave_eeg,'eeg');
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' EEG wavelet - average by scale'];
title(tt)

% --------------------------------------------------------------------
function varargout = unit_average_submenu_Callback(h, eventdata, handles, varargin)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_unit = get_unit_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Calculate averages and plot
wavelet_average(handles,wave_unit,'unit');
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' Unit wavelet - average by scale'];
title(tt)



% --------------------------------------------------------------------
% WAVELET AVERAGE subfunction
% --------------------------------------------------------------------
function wavelet_average(handles,wave,st)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Frequency bands
global SCALE
f = SCALE;
global WAVETIME
wavetime = WAVETIME;
fnd = find(f<1);
pwind1 = fnd(1);
fnd = find(f<3);
pwind2 = fnd(1);
fnd = find(f<6);
pwind3 = fnd(1);
fnd = find(f<20);
pwind4 = fnd(1);
fnd = find(f<50);
pwind5 = fnd(1);

% Creating "wavelet vectors"
wavepower = abs(wave).^2;
lenw = size(wave,2);
average = zeros(4,lenw);                        %average contains the mean along scale of five parts of wavelet power
average(1,:) = mean(wavepower(pwind5:pwind4-1,:));      %20 - 50 Hz
average(2,:) = mean(wavepower(pwind4:pwind3-1,:));     %6 - 20 Hz
average(3,:) = mean(wavepower(pwind3:pwind2-1,:));    %3 - 6 Hz
average(4,:) = mean(wavepower(pwind2:pwind1-1,:));   %1 - 3 Hz

% Plot
H = figure;
hold on
plot(average(1,:),'r');
plot(average(2,:),'b');
plot(average(3,:),'g');
plot(average(4,:),'m');
legend('20 - 50 Hz','6 - 20 Hz','3 - 6','1 - 3 Hz',2);
hold off
y_lim = ylim;
intlen =  length(average);
axis([0 intlen y_lim(1) y_lim(2)]);
b_rescaleaxis('X',wavetime)
setappdata(gca,'scalex',wavetime)
b_zoomset_for_wavelet

% Store wavelet in handles structure
nm = ['average_' st];
nm = upper(nm);
nmcoord = [nm '_coord'];
nmcoord = upper(nmcoord);
eval(['global ' nm]);
eval([nm ' =  average;']);
eval(['global ' nmcoord]);
eval([nmcoord ' =  [datinx1 datinx2];']);



% --------------------------------------------------------------------
% Average fft submenu callbacks
% --------------------------------------------------------------------
function eeg_average_fft_submenu_Callback(h, eventdata, handles, varargin)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_eeg = get_eeg_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Check if the appropriate wavelet average already exist
average_eeg = get_eeg_average(h,eventdata,handles,varargin,datinx1,datinx2);

% Calculate fft and plot
average_fft(handles,wave_eeg,average_eeg);
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' EEG average fft'];
title(tt)
setappdata(gcf,'title',tt);

% Set keypress function
set(gcf,'KeypressFcn','b_analyse4(''b_callfftswitch'',''eeg'')');

% --------------------------------------------------------------------
function unit_average_fft_submenu_Callback(h, eventdata, handles, varargin)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_unit = get_unit_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Check if the appropriate wavelet average already exist
average_unit = get_unit_average(h,eventdata,handles,varargin,datinx1,datinx2);

% Calculate fft and plot
average_fft(handles,wave_unit,average_unit);
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' Unit average fft'];
title(tt)
setappdata(gcf,'title',tt);

% Set keypress function
set(gcf,'KeypressFcn','b_analyse4(''b_callfftswitch'',''unit'')');



% --------------------------------------------------------------------
% AVERAGE FFT subfunction
% --------------------------------------------------------------------
function average_fft(handles,wave,average)

% Sampling rate
wavprop = handles.wavelet_properties;
samprate = wavprop.samprate;

% Calculate fft using Welch's method
lenwv = size(average,2);

global WHICHPERIODOGRAM
whichperiod = WHICHPERIODOGRAM;

H = figure;
switch whichperiod
case 1
    [p_1 f_1] = b_fft2(average(4,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it1 - 3 Hz}')
case 2
    [p_2 f_2] = b_fft2(average(3,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it3 - 6 Hz}')
case 3
    [p_3 f_3] = b_fft2(average(2,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it6 - 20 Hz}')
case 4
    [p_4 f_4] = b_fft2(average(1,:),samprate,lenwv,2);
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
    text(xind,0.5*ym,'{\it20 - 50 Hz}')
end

% Set keypress functions
global HANDLES
HANDLES = handles;
set(H,'KeyPressFcn','b_callfftswitch')



% --------------------------------------------------------------------
% Wavelet fft submenu callbacks
% --------------------------------------------------------------------
function varargout = eeg_wavelet_fft_submenu_Callback(h, eventdata, handles, varargin)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_eeg = get_eeg_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Import 'f'
global SCALE
f = SCALE;

% Sampling rate
wavprop = handles.wavelet_properties;

% Wavaletfft
b_waveletfft2(wave_eeg,f,wavprop.samprate);

% Title
axs = findobj(gcf,'Type','axes');
axes(axs(3))
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' EEG modulogram'];
title(tt)

% --------------------------------------------------------------------
function varargout = unit_wavelet_fft_submenu_Callback(h, eventdata, handles, varargin)

% Gain datinx
datinx1 = handles.datinx1;
datinx2 = handles.datinx2;

% Check if the appropriate wavelets already exist
wave_unit = get_unit_wavelet(h,eventdata,handles,varargin,datinx1,datinx2);

% Import 'f'
global SCALE
f = SCALE;

% Sampling rate
wavprop = handles.wavelet_properties;

% Wavaletfft
b_waveletfft2(wave_unit,f,wavprop.samprate);

% Title
axs = findobj(gcf,'Type','axes');
axes(axs(3))
tt = [handles.fname(1:3) ' ' handles.fname(5:6) ' ' num2str(datinx1)...
        ' ' num2str(datinx2) ' Unit modulogram'];
title(tt)



% --------------------------------------------------------------------
% GET WAVELET subfunctions
% --------------------------------------------------------------------
function wave_eeg = get_eeg_wavelet(h,eventdata,handles,varargin,datinx1,datinx2)

% Check if the appropriate wavelets already exist
global WAVE_EEG
if isempty(WAVE_EEG)
    eegpower_submenu_Callback(h,eventdata,handles,varargin)
    global WAVE_EEG
    wave_eeg = WAVE_EEG;
else
    global WAVE_EEG_COORD
    coord = WAVE_EEG_COORD;
    if ~isequal([datinx1 datinx2],coord)
        eegpower_submenu_Callback(h,eventdata,handles,varargin)
        global WAVE_EEG
        wave_eeg = WAVE_EEG;
    else
        global WAVE_EEG_PROP
        prop = WAVE_EEG_PROP;
        global WAVPROP_CURRENT
        newprop = WAVPROP_CURRENT;
        if ~isequal(prop,newprop) & ~isempty(newprop)
            eegpower_submenu_Callback(h,eventdata,handles,varargin)
            global WAVE_EEG
            wave_eeg = WAVE_EEG;
        else
            wave_eeg = WAVE_EEG;
        end
    end
end

% --------------------------------------------------------------------
function wave_unit = get_unit_wavelet(h,eventdata,handles,varargin,datinx1,datinx2)

% Check if the appropriate wavelets already exist
global WAVE_UNIT
if isempty(WAVE_UNIT)
    unitpower_submenu_Callback(h,eventdata,handles,varargin)
    global WAVE_UNIT
    wave_unit = WAVE_UNIT;
else
    global WAVE_UNIT_COORD
    coord = WAVE_UNIT_COORD;
    if ~isequal([datinx1 datinx2],coord)
        unitpower_submenu_Callback(h,eventdata,handles,varargin)
        global WAVE_UNIT
        wave_unit = WAVE_UNIT;
    else
        global WAVE_UNIT_PROP
        prop = WAVE_UNIT_PROP;
        global WAVPROP_CURRENT
        newprop = WAVPROP_CURRENT;
        if ~isequal(prop,newprop) & ~isempty(newprop)
            unitpower_submenu_Callback(h,eventdata,handles,varargin)
            global WAVE_UNIT
            wave_unit = WAVE_UNIT;
        else
            wave_unit = WAVE_UNIT;
        end
    end
end



% --------------------------------------------------------------------
% GET AVERAGE subfunctions
% --------------------------------------------------------------------
function average_eeg = get_eeg_average(h, eventdata, handles, varargin,datinx1,datinx2)

% Check if the appropriate wavelet average already exist
global AVERAGE_EEG
if isempty(AVERAGE_EEG)
    eeg_average_submenu_Callback(h,eventdata,handles,varargin)
    global AVERAGE_EEG
    average_eeg = AVERAGE_EEG;
else
    global AVERAGE_EEG_COORD
    coord = AVERAGE_EEG_COORD;
    if ~isequal([datinx1 datinx2],coord)
        eeg_average_submenu_Callback(h, eventdata, handles, varargin)
        global AVERAGE_EEG
        average_eeg = AVERAGE_EEG;
    else
        average_eeg = AVERAGE_EEG;
    end
end

% --------------------------------------------------------------------
function average_unit = get_unit_average(h, eventdata, handles, varargin,datinx1,datinx2)

% Check if the appropriate wavelet average already exist
global AVERAGE_UNIT
if isempty(AVERAGE_UNIT)
    unit_average_submenu_Callback(h, eventdata, handles, varargin)
    global AVERAGE_UNIT
    average_unit = AVERAGE_UNIT;
else
    global AVERAGE_UNIT_COORD
    coord = AVERAGE_UNIT_COORD;
    if ~isequal([datinx1 datinx2],coord)
        unit_average_submenu_Callback(h, eventdata, handles, varargin)
        global AVERAGE_UNIT
        average_unit = AVERAGE_UNIT;
    else
        average_unit = AVERAGE_UNIT;
    end
end



% --------------------------------------------------------------------
% Setup submenu callback
% --------------------------------------------------------------------
function varargout = setup_submenu_Callback(h, eventdata, handles, varargin)

% Launch WAVELET_SETTINGS GUI
b_wavelet_settings;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      S U B F U N C T I O N S                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------------------------------------------------------------
% Set slider subfunction
% --------------------------------------------------------------------
function set_Slider(h, eventdata, handles, varargin)

% Set slider
time = handles.time;
xxx = get(handles.axes_eeg,'XLim');
windowwidth = xxx(2) - xxx(1);
handles.windowwidth = windowwidth;
guidata(h,handles);
slidstep = windowwidth / (time(end) - windowwidth + eps - time(1));
slidstep_default = (time(end) - time(1)) / (time(end) - windowwidth + eps - time(1));
if slidstep == slidstep_default
    set(handles.slider1,'Enable','off');
else
    set(handles.slider1,'Value',xxx(1),'Min',time(1),'Max',time(end)-windowwidth+eps,'SliderStep',...
        [0.01,slidstep],'Enable','on');
end
drawnow;



% --------------------------------------------------------------------
% Set edit strings subfunction
% --------------------------------------------------------------------
function set_EditStrings(h, eventdata, handles, varargin)

% Set edit strings
x_lim = get(handles.axes_eeg,'XLim');
x_lim = round(x_lim*10) / 10;
set(handles.edit1,'String',num2str(x_lim(1)));
set(handles.edit2,'String',num2str(x_lim(2)));
handles = guidata(h);



% --------------------------------------------------------------------
% Set button down functions subfunction
% --------------------------------------------------------------------
function set_ButtonDownFcns(h, eventdata, handles, varargin)

% Set buttondown functions
bdf1 = 'b_analyse4(''axes_eeg_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
bdf2 = 'b_analyse4(''axes_unit_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
ach1 = get(handles.axes_eeg,'Children');
ach2 = get(handles.axes_unit,'Children');
set(handles.axes_eeg,'ButtonDownFcn',bdf1);
set(handles.axes_unit,'ButtonDownFcn',bdf2);
set(ach1,'ButtonDownFcn',bdf1);
set(ach2,'ButtonDownFcn',bdf2);



% --------------------------------------------------------------------
% Set datinx subfunction
% --------------------------------------------------------------------
function setdatinx(h,handles)

% Gain datinx
x_lim(1) = str2num(get(handles.edit1,'String'));
x_lim(2) = str2num(get(handles.edit2,'String'));
len = length(handles.time);
datinx1 = round(x_lim(1)*handles.sr+1);
datinx2 = round(x_lim(2)*handles.sr);
datinx2 = min(datinx2,len);

% Set datinx
handles.datinx1 = datinx1;
handles.datinx2 = datinx2;
guidata(h,handles)



% -------------------------------------------------------------------------
% Set focus subfunction
% -------------------------------------------------------------------------
function setfocus(h, handles, focus)

% Plot
y_eeg = get(handles.axes_eeg,'YLim');
y_unit = get(handles.axes_unit,'YLim');
time = handles.time;

axes(handles.axes_eeg);
axis([max(focus(1)-0.05,time(1)) min(focus(2)+0.05,time(end)) y_eeg(1) y_eeg(2)]);
b_liner(handles.plot1_handle,handles.eegs,handles.gui_pos,handles.la,handles.len,handles.sr);

axes(handles.axes_unit);
axis([max(focus(1)-0.05,time(1)) min(focus(2)+0.05,time(end)) y_eeg(1) y_eeg(2)]);
b_liner2(handles.plot2_handle,handles.units,handles.gui_pos,handles.la,handles.len,handles.sr);

drawnow;

% Set slider
eventdata = [];
set_Slider(h, eventdata, handles)
handles = guidata(h);

% Set edit strings
set_EditStrings(h, eventdata, handles)
handles = guidata(h);

% Set datinx
setdatinx(h,handles)
handles = guidata(h);