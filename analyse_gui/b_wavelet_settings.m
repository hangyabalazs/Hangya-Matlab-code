function varargout = b_wavelet_settings(varargin)
% B_WAVELET_SETTINGS Application M-file for b_wavelet_settings.fig
%    FIG = B_WAVELET_SETTINGS launch b_wavelet_settings GUI.
%    B_WAVELET_SETTINGS('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 31-Aug-2004 14:46:42

%WAVELET_SETTINGS   Dependent GUI of ANALYSE.
%   Wavelet settings can be modified through WAVELET_SETTINGS GUI. This GUI is
%   launched when user chooses 'Analysis' - 'Wavelet' - 'Settings' submenu in
%   ANALYSE GUI. The set values are remembered until next reset. See WAVELET for
%   a detailed help on wavelet parameters.
%
%   See also ANALYSE and WAVELET.

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
    
    % Load last setting
    global DATAPATH
    try
        ff = fullfile(DATAPATH,'analyse_gui\settings\wavprop.mat');
    end
    if b_isfilename(ff)
        load(ff)
        set(handles.samprate,'String',wavprop.samprate);
        set(handles.mif,'String',wavprop.mif);
        set(handles.mother,'String',wavprop.mother);
        set(handles.param,'String',wavprop.param);
        set(handles.dj,'String',wavprop.dj);
        set(handles.cfd,'Value',wavprop.cfd);
        set(handles.smj1,'Value',wavprop.smj1);
        if ~isempty(wavprop.j1)
            set(handles.j1,'Enable','on');
            set(handles.j1,'String',wavprop.j1);
        end
        set(handles.dst,'Value',wavprop.dst);
        set(handles.sms0,'Value',wavprop.sms0);
        if ~isempty(wavprop.s0)
            set(handles.s0,'Enable','on');
            set(handles.s0,'String',wavprop.s0);
        end
        set(handles.yes,'Value',wavprop.yes);
        set(handles.no,'Value',wavprop.no);
    else
        warndlg('Cannot load wavelet settings: permission denied or file does not exist. Setting values to default.',...
            'Loading failure');
        lasterr = [];
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
        if ~isempty(lasterr)
		    disp(lasterr);
            errordlg(lasterr,'Error','modal');
        end
	end

end



% --------------------------------------------------------------------
function varargout = samprate_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = mif_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = dj_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = s0_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = j1_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = mother_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
function varargout = param_Callback(h, eventdata, handles, varargin)



% --------------------------------------------------------------------
% DEFAULT button callback - reset values
% --------------------------------------------------------------------
function varargout = default_Callback(h, eventdata, handles, varargin)
set(handles.samprate,'String','400');
set(handles.mif,'String','0.5');
set(handles.mother,'String','morlet');
set(handles.param,'String','6');
set(handles.dj,'String','0.05');
set(handles.cfd,'Value',1);
set(handles.smj1,'Value',0);
set(handles.j1,'String','');
set(handles.j1,'Enable','off');
set(handles.dst,'Value',1);
set(handles.sms0,'Value',0);
set(handles.s0,'String','');
set(handles.s0,'Enable','off');
set(handles.yes,'Value',1);
set(handles.no,'Value',0);
% IN CASE YOU MODIFY DEFAULT SETTINGS, YOU HAVE TO MODIFY THEM IN ANALYSE.M
% LAUNCH PART AS WELL!



% --------------------------------------------------------------------
% CANCEL button callback - delete GUI
% --------------------------------------------------------------------
function varargout = cancel_Callback(h, eventdata, handles, varargin)
delete(handles.figure1)




% --------------------------------------------------------------------
% OK button callback - save settings and delete GUI
% --------------------------------------------------------------------
function varargout = ok_Callback(h, eventdata, handles, varargin)

% Gain GUI values and check if they have been set properly
wavprop.samprate = get(handles.samprate,'String');
wavprop.mif = get(handles.mif,'String');
wavprop.param = get(handles.param,'String');
wavprop.dj = get(handles.dj,'String');

fns = fieldnames(wavprop);
lfns = length(fns);
is = zeros(1,lfns);      % 'is': are there empty fields
iis = zeros(1,lfns);     % 'iis': are there non numeric values
for k = 1:lfns
    s = ['wavprop.' fns{k}];
    is(k) = eval(['isempty(' s ');']);
    eval([s ' = str2num(' s ');']);
    iis(k) = eval(['isempty(' s ');']);
end

wavprop.mother = get(handles.mother,'String');
if isempty(find(strcmp(wavprop.mother,{'morlet','paul','dog'})))
    warndlg('Invalid string for ''mother''. Set ''mother'' to ''morlet'', ''paul'' or ''dog''.','Warning');
    return
end

is(end+1) = isempty(wavprop.mother);

wavprop.cfd = get(handles.cfd,'Value');
wavprop.smj1 = get(handles.smj1,'Value');
wavprop.j1 = get(handles.j1,'String');

if wavprop.smj1
    is(end+1) = isempty(wavprop.j1);
    wavprop.j1 = str2num(wavprop.j1);
    iis(end+1) = isempty(wavprop.j1);
end
       
wavprop.dst = get(handles.dst,'Value');
wavprop.sms0 = get(handles.sms0,'Value');
wavprop.s0 = get(handles.s0,'String');

if wavprop.sms0
    is(end+1) = isempty(wavprop.s0);
    wavprop.s0 = str2num(wavprop.s0);
    iis(end+1) = isempty(wavprop.s0);
end

wavprop.yes = get(handles.yes,'Value');
wavprop.no = get(handles.no,'Value');

fis = find(is);
if ~isempty(fis)
    warndlg('Some of the parameters have not been set.','Warning');
    return
end
fiis = find(iis);
if ~isempty(fiis)
    warndlg('Numeric parameter has been set to non numeric value.','Warning');
    return
end

% Save settings and delete GUI
global DATAPATH
ff = fullfile(DATAPATH,'analyse_gui\settings\wavprop.mat');
str = ['save ' ff ' wavprop'];
eval(str)
global WAVPROP_CURRENT
WAVPROP_CURRENT = wavprop;
delete(handles.figure1)



% --------------------------------------------------------------------
% 'CALCULATE FROM DJ' radiobutton callback
% --------------------------------------------------------------------
function varargout = cfd_Callback(h, eventdata, handles, varargin)
set(handles.cfd,'Value',1)
set(handles.smj1,'Value',0)
set(handles.j1,'String','')
set(handles.j1,'Enable','off')



% --------------------------------------------------------------------
% 'SET J1 MANUALLY' radiobutton callback
% --------------------------------------------------------------------
function varargout = smj1_Callback(h, eventdata, handles, varargin)
set(handles.smj1,'Value',1)
set(handles.cfd,'Value',0)
set(handles.j1,'Enable','on')



% --------------------------------------------------------------------
% 'DOUBLE SAMPLING TIME' radiobutton callback
% --------------------------------------------------------------------
function varargout = dst_Callback(h, eventdata, handles, varargin)
set(handles.dst,'Value',1)
set(handles.sms0,'Value',0)
set(handles.s0,'String','')
set(handles.s0,'Enable','off')



% --------------------------------------------------------------------
% 'SET S0 MANUALLY' radiobutton callback
% --------------------------------------------------------------------
function varargout = sms0_Callback(h, eventdata, handles, varargin)
set(handles.sms0,'Value',1)
set(handles.dst,'Value',0)
set(handles.s0,'Enable','on')



% --------------------------------------------------------------------
% 'YES' radiobutton callback
% --------------------------------------------------------------------
function varargout = yes_Callback(h, eventdata, handles, varargin)
set(handles.yes,'Value',1)
set(handles.no,'Value',0)



% --------------------------------------------------------------------
% 'NO' radiobutton callback
% --------------------------------------------------------------------
function varargout = no_Callback(h, eventdata, handles, varargin)
set(handles.no,'Value',1)
set(handles.yes,'Value',0)