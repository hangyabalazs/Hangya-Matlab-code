function somplotHindex
%SOMPLOTHINDEX   Plots H-index for interneuron paper.
%
%   See also SOMISSTIM3.

% Load Cellbase
loadcb

% Select upon clustering criteria
ID = getvalue('ID',CELLIDLIST);
Lratio = getvalue('Lratio',CELLIDLIST);
% global cells
cells = CELLIDLIST((ID>18|isnan(ID))&(Lratio<0.15|isnan(Lratio)));

% Plot H-index
Hindex = getvalue('H_index',cells);
figure
hist(Hindex,1000)
P = findobj(gca,'Type','patch');
set(P,'FaceColor',[0.75 0 0],'EdgeColor',[0.75 0 0])
set(gca,'Box','off','LineWidth',2,'FontSize',16,'FontWeight','bold',...
    'XScale','log','YTick',[0 50 100],'XLim',[0.0007 1],...
    'XMinorTick','off','TickDir','out')

cd('c:\Balazs\_analysis\SOM\KL\somisstim3\solution_for_inhibition\')
close all
s1 = find(Hindex<0.01);
for k = 4:length(s1)
    disp(s1(k))
    disp(cells{s1(k)})
    [p_value Idiff p_value2 Idiff2 p_value3 Idiff3 p_value4 Idiff4] = ...
        somisstim3(cells{s1(k)});
%     saveas(H,[cells{s1(k)} '_slsi.fig'])
%     disp(p_value4)
%     close all
end

s2 = find(Hindex>=0.01&Hindex<0.05);
for k = 1:length(s2)
    disp(s2(k))
    disp(cells{s2(k)})
    [p_value Idiff p_value2 Idiff2 p_value3 Idiff3 p_value4 Idiff4] = ...
        somisstim3(cells{s2(k)});
%     saveas(H,[cells{s1(k)} '_slsi.fig'])
%     disp(p_value4)
%     close all
end

keyboard
    