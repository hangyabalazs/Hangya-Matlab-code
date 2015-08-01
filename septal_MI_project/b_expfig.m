function b_expfig
%EXPFIG   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the figure in a new figure window (exports only spike train,
%           return plot or variance plot)
%
%   See also EXPORTFIG, CALLEXPORTFIG and CALLEXPFIG.

% Whitch plot to export
c = input('Spike Train: 1, Return Plot: 2, Variance Plot: 3 ');

% Get handles structure
global HANDLES
handles = HANDLES;
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
d = MaxClusterNumber;

index_selected2 = get(handles.listbox2,'Value');
index_selected1 = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
fln = file_list{index_selected1};
i_first = str2num(fln(11:16));
i_second = str2num(fln(18:23));

% Export
figure
switch c
case 1
    z = zeros(1,length(Time));
    Vdisc_new = Vdisc - (i_first - 320000);
    z(Vdisc_new) = 1;
    z(Vdisc_new+1) = -1;
    dec = index_selected2 + 1;
%     plot2_handle = plot(Time,z,'Color',[ 0.631 0.941 1.000 ]);
    plot2_handle = plot(Time,z,'r');
    axis([Time(1) Time(end) -1.5 1.5])
    hold on;
    for j = 1:size(Burst{dec},2)
        ind1 = Vdisc_new(Burst{dec}(1,j));
        ind2 = Vdisc_new(Burst{dec}(2,j));
        plot(Time(ind1-1:ind2+1),z(ind1-1:ind2+1),'b')
%         rajz = [Time(Vdisc_new(Burst{dec}(1,j))) Time(Vdisc_new(Burst{dec}(2,j)));0.2 0.7];
%         line_handle = line(rajz(1,:),rajz(2,:),'Color','k');
%         set(line_handle,'ButtonDownFcn',bdf2)
    end
    hold off
case 2
    plot(ReturnPlotXData,ReturnPlotYData,'.');
case 3
    hold on
    plot([1:d],IntraBurstIvVar,'g');
    plot([1:d],ExtraBurstIvVar,'b');
    plot([1:d],InterBurstIvVar,'c');
    plot([1:d],FirstSpikeVar,'m');
    plot([1:d],AllFirstSpikeVar,'k');
    legend('intraburstivvar','extraburstivvar','interburstivvar','firstspikevar','allfirstspikevar',2);
    xlim([1 d]);
    x_lim = xlim;
    y_lim = ylim;
    kk = (y_lim(2) - y_lim(1)) / 3;
    hold off
end