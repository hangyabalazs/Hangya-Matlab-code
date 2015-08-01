function b_zoom_for_wavelet_linkaxes
%ZOOM_FOR_WAVELET_LINKAXES  Zoom in and out on a 2-D plot.
%   ZOOM_FOR_WAVELET_LINKAXES acts similar to ZOOM. You can modify each
%   step individualy and assign ZOOM_FOR_WAVELET_LINKAXES as figure
%   ButtonDownFcn with ZOOMSET_FOR_WAVELET.
%
%   ZOOM_FOR_WAVELET_LINKAXES 'rescales' figure after each zooming. Note,
%   that you have to define the scale vectors as application data, so
%   include the following lines in your program:
%
%       setappdata(gca,'scalex',wavetime)
%       setappdata(gca,'scaley',f)
%       b_zoomset_for_wavelet
%
%   where 'wavetime' and 'f' are the scale vectors.
%
%   ZOOM_FOR_WAVELET_LINKAXES performs linked x-axis zoom in another
%   subplot of the figure. The handle of the other subplot should be given
%   as application data as follows:
%
%       setappdata(gca,'subplothandle',S2)
%
%   where 'S2' is the handle of the subplot to be linked.
%
%   See also ZOOM, ZOOM2, ZOOMSET_FOR_WAVELET, ZOOMSET, RETIMEFIG and RESCALEFIG.

% Get handles
fig = gcf;
ax = gca;
plt = findobj(ax,'Type','line');
ptch = findobj(ax,'Type','patch');
im = findobj(ax,'Type','image');

% Get application data
x_inidata = getappdata(ax,'x_data');
y_inidata = getappdata(ax,'y_data');
if ~iscell(x_inidata)
    tempcell{1} = x_inidata;
    x_inidata = tempcell;
    tempcell{1} = y_inidata;
    y_inidata = tempcell;
end
x_inilim = getappdata(ax,'x_lim');
y_inilim = getappdata(ax,'y_lim');

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
    set(ax,'XTickLabelMode','auto')     % set TickLabelMode
    set(ax,'YTickLabelMode','auto')
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
        for n = 1:length(x_inidata)
            xxa = find(x_inidata{n}<=xn(1));
            xxb = find(x_inidata{n}>=xn(2));
            if ~isempty(xxa)
                xxn(1) = xxa(end);
            else
                xxn(1) = 1;
            end
            if ~isempty(xxb)
                xxn(2) = xxb(1);
            else
                xxn(2) = length(x_inidata{n});
            end
            set(plt(n),'YData',y_inidata{n}(xxn(1):xxn(2)),'XData',x_inidata{n}(xxn(1):xxn(2)));
        end
    end
    
case 'open'     % reset original plot
    set(ax,'XTickLabelMode','auto')     % set TickLabelMode
    set(ax,'YTickLabelMode','auto')
    if ~isempty(plt)
        for n = 1:length(x_inidata)
            set(plt(n),'YData',y_inidata{n},'XData',x_inidata{n});
        end
    elseif ~isempty(ptch) | ~isempty(im)
        set(ax,'YLim',y_inilim,'XLim',x_inilim);
    end
    axis([x_inilim(1) x_inilim(2) y_inilim(1) y_inilim(2)]);
    
otherwise   % zoom out
    set(ax,'XTickLabelMode','auto')     % set TickLabelMode
    set(ax,'YTickLabelMode','auto')
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
        for n = 1:length(x_inidata)
            xxa = find(x_inidata{n}<=xn(1));
            xxb = find(x_inidata{n}>=xn(2));
            if ~isempty(xxa)
                xxn(1) = xxa(end);
            else
                xxn(1) = 1;
            end
            if ~isempty(xxb)
                xxn(2) = xxb(1);
            else
                xxn(2) = length(x_inidata{n});
            end
            set(plt(n),'YData',y_inidata{n}(xxn(1):xxn(2)),'XData',x_inidata{n}(xxn(1):xxn(2)));
        end
    end
end

% Rescale
scalex = getappdata(ax,'scalex');
scaley = getappdata(ax,'scaley');

if ~isempty(scalex)
    b_rescaleaxis('X',scalex)
end
if ~isempty(scaley)
    b_rescaleaxis('Y',scaley)
end

drawnow;

% Link axes
xxn = xlim(ax);
subplothandle = getappdata(ax,'subplothandle');
for k = 1:length(subplothandle)
    xlim(subplothandle(k),xxn)
end