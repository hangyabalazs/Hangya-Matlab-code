function b_greyplot(s)
%GREYPLOT   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       y - plots wavelet average magnitudes from 1 to 3 hz against time on the spike train
%       x - plots wavelet average magnitudes from 3 to 6 hz against time on the spike train
%       c - plots wavelet average magnitudes from 6 to 20 hz against time on the spike train
%       v - plots wavelet average magnitudes from 20 to 50 hz against time on the spike train
%
%   See also CALLGREYPLOT.

% Get handles structure
global HANDLES
handles = HANDLES;
file_list = get(handles.listbox1,'String');
index_selected1 = get(handles.listbox1,'Value');
global DATAPATH
ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave\',file_list{index_selected1});
load(ff)

sw = size(Wavevec,2);
tm = sw * 100;
ttm = linspace(0,tm,sw);

% Plot
subplot(handles.subplot2_handle);
hold on
switch s
case 'r'
    if isempty(handles.greyplot1_handle)
        wm = Wavevec(1,:) - mean(Wavevec(1,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot1_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot1_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot1_handle,'Visible','off')
        else
            set(handles.greyplot1_handle,'Visible','on')
            if ~isempty(handles.greyplot2_handle)
                set(handles.greyplot2_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot3_handle)
                set(handles.greyplot3_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot4_handle)
                set(handles.greyplot4_handle,'Visible','off')
            end
        end
    end
case 'b'
    if isempty(handles.greyplot2_handle)
        wm = Wavevec(2,:) - mean(Wavevec(2,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot2_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot2_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot2_handle,'Visible','off')
        else
            set(handles.greyplot2_handle,'Visible','on')
            if ~isempty(handles.greyplot1_handle)
                set(handles.greyplot1_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot3_handle)
                set(handles.greyplot3_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot4_handle)
                set(handles.greyplot4_handle,'Visible','off')
            end
        end
    end
case 'g'
    if isempty(handles.greyplot3_handle)
        wm = Wavevec(3,:) - mean(Wavevec(3,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot3_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot4_handle)
            set(handles.greyplot4_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot3_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot3_handle,'Visible','off')
        else
            set(handles.greyplot3_handle,'Visible','on')
            if ~isempty(handles.greyplot1_handle)
                set(handles.greyplot1_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot2_handle)
                set(handles.greyplot2_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot4_handle)
                set(handles.greyplot4_handle,'Visible','off')
            end
        end
    end
case 'm'
    if isempty(handles.greyplot4_handle)
        wm = Wavevec(4,:) - mean(Wavevec(4,:));
        wwm = wm / max(wm) * 1.5;
        handles.greyplot4_handle = plot(ttm,wwm,'k');
        if ~isempty(handles.greyplot1_handle)
            set(handles.greyplot1_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot2_handle)
            set(handles.greyplot2_handle,'Visible','off')
        end
        if ~isempty(handles.greyplot3_handle)
            set(handles.greyplot3_handle,'Visible','off')
        end
    else
        vis = get(handles.greyplot4_handle,'Visible');
        if strcmp(vis,'on')
            set(handles.greyplot4_handle,'Visible','off')
        else
            set(handles.greyplot4_handle,'Visible','on')
            if ~isempty(handles.greyplot1_handle)
                set(handles.greyplot1_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot2_handle)
                set(handles.greyplot2_handle,'Visible','off')
            end
            if ~isempty(handles.greyplot3_handle)
                set(handles.greyplot3_handle,'Visible','off')
            end
        end
    end
end
HANDLES = handles;
global HANDLES
hold off