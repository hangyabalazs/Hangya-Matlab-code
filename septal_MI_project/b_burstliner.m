function b_burstliner
%BURSTLINER Keypress function for ICA_GUI2B.
%   Pressing 'l' in the figure window of 'Iterative Cluster Analysis' graphical user
%   interface, lines connecting the actual bursts become visible. Pressing 'l' again,
%   the lines disappear.
%
%   See also CALLGREYPLOT.

% Loadig ICA GUI informations
global HANDLES
handles  = HANDLES;
index_selected1 = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
index_selected2 = get(handles.listbox2,'Value');
dec = index_selected2 + 1;

% Loadig cluster analysis data
global DATAPATH
cl = file_list{index_selected1}(12:31);
pth = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly'];
fln = ['THETA_ICA_' cl];
ff = fullfile(pth,fln);
load(ff)
i_first = str2num(fln(18:23));
i_second = str2num(fln(25:30));

% Drawing the lines or switching their visibility
hold on;
sbd = size(Burst{dec},2);
Vdisc_new = Vdisc - (i_first - 320000);
if isempty(handles.line_handles)
    line_handles = zeros(1,sbd);
    for j = 1:sbd
        rajz = [Time(Vdisc_new(Burst{dec}(1,j))) Time(Vdisc_new(Burst{dec}(2,j)));0.2 0.7];
        line_handles(j) = line(rajz(1,:),rajz(2,:),'Color','k');
    end
    handles.line_handles = line_handles;
else
    vis = get(handles.line_handles,'Visible');
    if strcmp(vis{1},'on')
        set(handles.line_handles,'Visible','off')
    else
        set(handles.line_handles,'Visible','on')
    end
end
hold off

% Modifying ICA GUI
HANDLES = handles;
global HANDLES