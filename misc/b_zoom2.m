function zoom2
%ZOOM2  Zoom in and out on a 2-D plot.
%   ZOOM2 acts similar to ZOOM. You can modify each step individualy
%   and assign ZOOM2 as figure ButtonDownFcn with ZOOMSET.
%
%   See also ZOOM and ZOOMSET.

% Get handles
fig = gcf;
ax = gca;
plt = findobj(fig,'Type','line');
ptch = findobj(fig,'Type','patch');

% Get application data
x_inidata = getappdata(fig,'x_data');
y_inidata = getappdata(fig,'y_data');
x_inilim = getappdata(fig,'x_lim');
y_inilim = getappdata(fig,'y_lim');

% Zoom
seltyp = get(fig,'SelectionType');
switch seltyp
case 'normal'   % zoom in
    point1 = get(ax,'CurrentPoint'); % button down detected
    units = get(fig,'units');
    set(fig,'units','pixels')
    rbbox([get(fig,'currentpoint') 0 0],get(fig,'currentpoint'),fig);                   % return figure units
    set(fig,'units',units)
    point2 = get(ax,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    if isequal(point1,point2),  % set new axis limits
        xx = get(ax,'Xlim');
        yy = get(ax,'Ylim');
        xx2 = (abs(xx(2) - xx(1))) / 4;
        yy2 = (abs(yy(2) - yy(1))) / 4;
        xx3(1) = point1(1) - xx2;
        xx3(2) = point1(1) + xx2;
        yy3(1) = point1(2) - yy2;
        yy3(2) = point1(2) + yy2;
        if x_inilim(1) < xx3(1) & x_inilim(2) > xx3(2),
            set(ax,'Xlim',xx3);
        elseif x_inilim(1) > xx3(1),
            xx_new(1) = x_inilim(1);
            xx_new(2) = x_inilim(1) + (2 * xx2);
            set(ax,'Xlim',xx_new);
        elseif x_inilim(2) < xx3(2),
            xx_new(1) = x_inilim(2) - (2 * xx2);
            xx_new(2) = x_inilim(2);
            set(ax,'Xlim',xx_new);
        end
        if y_inilim(1) < yy3(1) & y_inilim(2) > yy3(2),
            set(ax,'Ylim',yy3);
        elseif y_inilim(1) > yy3(1),
            yy_new(1) = y_inilim(1);
            yy_new(2) = y_inilim(1) + (2 * yy2);
            set(ax,'Ylim',yy_new);
        elseif y_inilim(2) < yy3(2),
            yy_new(1) = y_inilim(2) - (2 * yy2);
            yy_new(2) = y_inilim(2);
            set(ax,'Ylim',yy_new);
        end
    else
        m1 = min([point1(1) point2(1)]);
        m2 = max([point1(1) point2(1)]);
        if isequal(m1,m2)
            m2 = m2 + 0.0001;
        end
        m3 = min([point1(2) point2(2)]);
        m4 = max([point1(2) point2(2)]);
        if isequal(m3,m4)
            m4 = m4 + 0.0001;
        end
        axis([m1 m2 m3 m4]);
    end

    if ~isempty(plt)
        xn = get(ax,'XLim');    % set new plot data
        xxa = find(x_inidata<=xn(1));
        xxb = find(x_inidata>=xn(2));
        if ~isempty(xxa)
            xxn(1) = xxa(end);
        else
            xxn(1) = 1;
        end
        if ~isempty(xxb)
            xxn(2) = xxb(1);
        else
            xxn(2) = length(x_inidata);
        end
        set(plt,'YData',y_inidata(xxn(1):xxn(2)),'XData',x_inidata(xxn(1):xxn(2)));
    end
    
case 'open'     % reset original plot
    if ~isempty(plt)
        set(plt,'YData',y_inidata,'XData',x_inidata);
    elseif ~isempty(ptch)
        set(ax,'YLim',y_inilim,'XLim',x_inilim);
    end
    axis([x_inilim(1) x_inilim(2) y_inilim(1) y_inilim(2)]);
    
otherwise   % zoom out
    point = get(ax,'CurrentPoint'); % button down detected
    point = point(1,1:2);
    xx = get(ax,'Xlim');    % set new axis limits
    yy = get(ax,'Ylim');
    xx2 = abs(xx(2) - xx(1));
    yy2 = abs(yy(2) - yy(1));
    if xx2 > (x_inilim(2) - x_inilim(1)) / 2,
        xx2 = (x_inilim(2) - x_inilim(1)) / 2;
    end
    if yy2 > (abs(y_inilim(2) - y_inilim(1))) / 2,
        yy2 = (abs(y_inilim(2) - y_inilim(1))) / 2;
    end
    xx3(1) = point(1) - xx2;
    xx3(2) = point(1) + xx2;
    yy3(1) = point(2) - yy2;
    yy3(2) = point(2) + yy2;
    if x_inilim(1) < xx3(1) & x_inilim(2) > xx3(2),
        set(ax,'Xlim',xx3);
    elseif x_inilim(1) > xx3(1),
        xx_new(1) = x_inilim(1);
        xx_new(2) = x_inilim(1) + (2 * xx2);
        set(ax,'Xlim',xx_new);
        set(ax,'Xlim',xx_new);
    elseif x_inilim(2) < xx3(2),
        xx_new(1) = x_inilim(2) - (2 * xx2);
        xx_new(2) = x_inilim(2);
        set(ax,'Xlim',xx_new);
    end
    if y_inilim(1) < yy3(1) & y_inilim(2) > yy3(2),
        set(ax,'Ylim',yy3);
    elseif y_inilim(1) > yy3(1),
        yy_new(1) = y_inilim(1);
        yy_new(2) = y_inilim(1) + (2 * yy2);
        set(ax,'Ylim',yy_new);
    elseif y_inilim(2) < yy3(2),
        yy_new(1) = y_inilim(2) - (2 * yy2);
        yy_new(2) = y_inilim(2);
        set(ax,'Ylim',yy_new);
    end
    
    if ~isempty(plt)
        xn = get(ax,'XLim');    % set new plot data
        xxa = find(x_inidata<=xn(1));
        xxb = find(x_inidata>=xn(2));
        if ~isempty(xxa)
            xxn(1) = xxa(end);
        else
            xxn(1) = 1;
        end
        if ~isempty(xxb)
            xxn(2) = xxb(1);
        else
            xxn(2) = length(x_inidata);
        end
        set(plt,'YData',y_inidata(xxn(1):xxn(2)),'XData',x_inidata(xxn(1):xxn(2)));
    end
    
end

drawnow;