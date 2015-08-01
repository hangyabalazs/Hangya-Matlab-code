function zoomset
%ZOOMSET    Assigns ZOOM2 as figure ButtonDownFcn.
%
%   See also ZOOM2.

% Get handles
fig = gcf;
ax = gca;
plt = findobj(fig,'Type','line');
ptch = findobj(fig,'Type','patch');

% Set application data
if ~isempty(plt)
    x_data = get(plt,'XData');
    y_data = get(plt,'YData');
    setappdata(fig,'x_data',x_data)
    setappdata(fig,'y_data',y_data)
end
x_lim = get(ax,'XLim');
y_lim = get(ax,'YLim');
setappdata(fig,'x_lim',x_lim)
setappdata(fig,'y_lim',y_lim)

% Set ButtonDownFcn
set(fig,'ButtonDownFcn','zoom4')
set(ax,'ButtonDownFcn','zoom4')
if ~isempty(plt)
    set(plt,'ButtonDownFcn','zoom4')
elseif ~isempty(ptch)
    set(ptch,'ButtonDownFcn','zoom4')
end