function b_fftswitch(s)
%FFTSWITCH   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       q - switches up one wavelet average magnitude fft plot
%       w - switches down one wavelet average magnitude fft plot
%
%   See also CALLFFTSWITCH.

% Switch
global WHICHPERIODOGRAM
whichperiod = WHICHPERIODOGRAM;
switch s
case 'u'
    whichperiod = whichperiod + 1;
    whichperiod = mod(whichperiod,4);
    if whichperiod == 0
        whichperiod = 4;
    end
case 'd'
    whichperiod = whichperiod - 1;
    whichperiod = mod(whichperiod,4);
    if whichperiod == 0
        whichperiod = 4;
    end
end
WHICHPERIODOGRAM = whichperiod;
global WHICHPERIODOGRAM

% Get handles structure
global HANDLES
handles = HANDLES;
set(handles.figure1,'HandleVisibility','on')
axes(handles.axes1)

% Load wavelet average magnitude vectors
index_selected1 = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
global DATAPATH
ff = fullfile(DATAPATH,'ICA\ica_gui2b\ica_wave\',file_list{index_selected1});
load(ff)
lenwv = size(Wavevec,2);

% Plot fft
pp1 = log2(lenwv);
pp2 = floor(pp1);
pp3 = 2^pp2;
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
% xtl = get(handles.axes1,'XTickLabel');
% xtl2 = str2num(xtl);
% xtl3 = xtl2 * 100 / pp3;
% xtl4 = num2str(xtl3);
% set(handles.axes1,'XTickLabel',xtl4)
set(handles.figure1,'HandleVisibility','callback')