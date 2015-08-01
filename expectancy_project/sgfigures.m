function sgfigures
%SGFIGURES   Modify figure style.
%   SGFIGURES modifies figure style of expectancy result figures for
%   publication.

% Directory
global DATAPATH
inpdir = [DATAPATH 'Expectancy\FzCzPz_SE\rose\'];
mm = pwd;
cd(inpdir)

% Modify rose diagrams
fl = dir(inpdir);
for k = 3:length(fl)
    open(fl(k).name)
    rosemod
    fn = [fl(k).name(1:end-4) '_mod.fig'];
    saveas(gcf,fn)
end
cd(mm)

% -------------------------------------------------------------------------
function rosemod

ach = allchild(gca);
tx = findobj(ach,'type','text');
set(tx,'FontWeight','bold','FontSize',16)
ln = findobj(ach,'type','line');
set(ln,'Color','black')

stx = get(tx,'String');
for k = 1:length(stx)
    switch stx{k}
        case {'30' '60' '120' '150' '210' '240' '300' '330'}
            delete(tx(k))
        case {'Cz 10%' 'Cz 37%' 'Cz 64%' 'Cz 91%' 'Fz 10%' 'Fz 37%' 'Fz 64%' 'Fz 91%' 'Pz 10%' 'Pz 37%' 'Pz 64%' 'Pz 91%'}
            delete(tx(k))
        case '  50'
            set(tx(k),'Position',[30.3519 6.7449 3047.9665])
        case '  100'
            set(tx(k),'Position',[70.8211 31.0264 3047.9665])
        case '  150'
            set(tx(k),'Position',[115.3372 53.9589 3047.9665])
        case '  200'
            set(tx(k),'Position',[162.5513 79.5894 3047.9665])
    end
end

set(ln(1),'LineWidth',2)
set(ln(2),'LineStyle','-','LineWidth',1)
delete(ln(3))
delete(ln(4))
set(ln(5),'LineStyle','-','LineWidth',1)
delete(ln(6))
delete(ln(7))
set(ln(9),'LineStyle','-','LineWidth',1)
set(ln(10),'LineStyle','-','LineWidth',1)
set(ln(11),'LineStyle','-','LineWidth',1)