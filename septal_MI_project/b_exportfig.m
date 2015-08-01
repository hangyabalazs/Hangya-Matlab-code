function b_exportfig(w)
%EXPORTFIG   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the upper figure in a new figure window
%       r - exports the lower figure in a new figure window
%       t - exports the upper accessoric subplot in a new figure window
%       z - exports the lower accessoric subplot in a new figure window
%
%   See also CALLEXPORTFIG.

% Get handles structure
global HANDLES
handles = HANDLES;
set(handles.figure1,'HandleVisibility','on')

% Export figure
switch w
case 'e'
    currentaxes = handles.axes1;
case 'r'
    currentaxes = handles.axes2;
case 't'
    currentaxes = handles.subplot1_handle;
case 'z'
    currentaxes = handles.subplot2_handle;
end
h = figure;
copyobj(currentaxes,h)
set(gca,'units','normalized')
set(gca,'position',[0.1300    0.1100    0.7750    0.8150])
set(handles.figure1,'HandleVisibility','callback')