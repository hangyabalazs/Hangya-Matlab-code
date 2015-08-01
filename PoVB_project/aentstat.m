% Collection of entropy statistics between bas and bic.

load('F:\balazs\_analysis\Andi\Ketxyl\EntrySTat\AUJ24_bas9_egesz_ENTRY.mat')
rUxy_bas=rUxy;
rUyx_bas=rUyx;
load('F:\balazs\_analysis\Andi\Ketxyl\EntrySTat\AUJ24_bic3_180-vegig_ENTRY.mat')
rUyx_bic=rUyx;
rUxy_bic=rUxy;

boxplot([rUxy_bas rUxy_bic rUyx_bas rUyx_bic],[zeros(size(rUxy_bas)) ones(size(rUxy_bic)) 2*ones(size(rUyx_bas)) 3*ones(size(rUyx_bic))],'labels',[{'bas. EEG->unit'} {'bic. EEG->unit'} {'bas. unit->EEG'} {'bic. unit->EEG'}])
[Wp_eu,Wh_eu] = b_ranksum2(rUxy_bic,rUxy_bas,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 4;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
[Wp_eu,Wh_eu] = b_ranksum2(rUyx_bic,rUyx_bas,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * 3 / 4;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')

saveas(gcf,'AUJ24_bas9-bic3_BOXuu')

load('F:\balazs\_analysis\Andi\Ketxyl\EntrySTat\AUJ24_bas9_egesz_ENTRY.mat')
rIxy_bas=rIxy;
load('F:\balazs\_analysis\Andi\Ketxyl\EntrySTat\AUJ24_bic3_180-vegig_ENTRY.mat')
rIxy_bic=rIxy;

boxplot([rIxy_bas rIxy_bic],[zeros(size(rIxy_bas)) ones(size(rIxy_bic))],'labels',[{'bas. MI'} {'bic. MI'}])
[Wp_eu,Wh_eu] = b_ranksum2(rIxy_bic,rIxy_bas,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')

saveas(gcf,'AUJ24_bas9-bic3_BOXii')