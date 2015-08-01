%% pair report - load

int = 17;
pyr = 20;
exp = 'hux044-day04-tr1-base-04';
type = 'positive';

global DATAPATH
inpdir = [DATAPATH 'Czurko\discriminated2\pos\' exp '\'];
cd(inpdir)
load('pfs.mat')
load('placemaps.mat')
load('compl_index.mat')
load('placecell_index.mat')
int2 = find(pci==int);
pyr2 = find(pci==pyr);
open(['autocorr' num2str(int2) '.fig'])
A_autoint = gca;
open(['autocorr' num2str(pyr2) '.fig'])
A_autopyr = gca;
open(['crosscorr_' num2str(min(pyr2,int2)) '_' num2str(max(pyr2,int2)) '.fig'])
A_cross = gca;
open(['normcrosscorr_' num2str(min(pyr2,int2)) '_' num2str(max(pyr2,int2)) '.fig'])
A_normcross = gca;

intmap = rhst{int};
iintmap = interp2_nonnan(intmap,5);
pyrmap = rhst{pyr};
ipyrmap = interp2_nonnan(pyrmap,5);
comp1 = C(int2,pyr2);
comp2 = (comp1 - C(int2,int2)) / C(pyr2,pyr2);
spcorr = Rmod(int2,pyr2);
if int2 > pyr2
    ip = 'P -> I';
else
    ip = 'I -> P';
end

%% generate report

Hout = figure;
subplot(3,3,1)
pcolor(nanpad(iintmap,1))
shading flat
cl = get(gca,'CLim');
set(gca,'CLim',[0 cl(2)])
colorbar
axis off

subplot(3,3,2)
pcolor(nanpad(ipyrmap,1))
shading flat
cl = get(gca,'CLim');
set(gca,'CLim',[0 cl(2)])
colorbar
axis off

S = subplot(3,3,4);
ach1 = findobj(allchild(A_autoint),'Type','line');
ach2 = findobj(allchild(A_autoint),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0])
copyobj(ach1,S)
copyobj(ach2,S)
set(gca,'XLim',[-200 200],'YTick',[],'YTickLabel',[])
title('int. autocorr.')

S = subplot(3,3,5);
ach1 = findobj(allchild(A_autopyr),'Type','line');
ach2 = findobj(allchild(A_autopyr),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0])
copyobj(ach1,S)
copyobj(ach2,S)
set(gca,'XLim',[-200 200],'YTick',[],'YTickLabel',[])
title('pyr. autocorr.')

S = subplot(3,3,7);
ach1 = findobj(allchild(A_cross),'Type','line');
ach2 = findobj(allchild(A_cross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0])
copyobj(ach1,S)
copyobj(ach2,S)
set(gca,'XLim',[-50 50],'YTick',[],'YTickLabel',[])
title([ip ' crosscorr.'])

S = subplot(3,3,8);
ach1 = findobj(allchild(A_normcross),'Type','line');
ach2 = findobj(allchild(A_normcross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0])
copyobj(ach1,S)
copyobj(ach2,S)
set(gca,'XLim',[-50 50],'YTick',[],'YTickLabel',[])
title('norm. crosscorr.')

set(gcf,'Unit','normalized')
str = {exp type ['R: ' num2str(spcorr)] ['C: ' num2str(comp2)]};

uicontrol('Style','text','Unit','normalized','Position',...
    [0.65 0.6 0.3 0.3],'FontSize',12,'HorizontalAlignment',...
    'left','String',str,'BackgroundColor',[1 1 1]);

%% Export

% set(gcf,'Unit','points','PaperUnits','points','PaperPositionMode','manual')
% pos = get(gcf,'Position');
% set(gcf,'Position',[pos(1) pos(2) pos(3)/pos(4)*1000 1000])
% set(gcf,'PaperPosition',[20 20 pos(3) pos(4)])
fh = ['-f' num2str(Hout)];
fname = [exp '_' num2str(int) '_' num2str(pyr)];
print(fh,'-dpdf',fname)
close all