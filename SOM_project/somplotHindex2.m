function somplotHindex2
%SOMPLOTHINDEX2   Plots H-index for interneuron paper.
%
%   See also SOMISSTIM3.

% Load Cellbase
loadcb

% Select
validinx = getvalue('valid',CELLIDLIST);
valid_neurons = CELLIDLIST(find(validinx));

% Plot H-index
Hindex = getvalue('H_index',valid_neurons);
figure
hist(Hindex,1000)
P = findobj(gca,'Type','patch');
set(P,'FaceColor',[0.75 0 0],'EdgeColor',[0.75 0 0])
set(gca,'Box','off','LineWidth',2,'FontSize',16,'FontWeight','bold',...
    'XScale','log','YTick',[0 100 200 300 400],'XLim',[0.0007 1],...
    'XMinorTick','off','TickDir','out')
keyboard