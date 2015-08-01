%% cells in behavior

selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells

%% ID, L-ratio

ID_NT = getvalue('ID_PC',[NT pChAT]);
ID_ChAT = getvalue('ID_PC',ChAT);

Lr_NT = getvalue('Lr_PC',[NT pChAT]);
Lr_ChAT = getvalue('Lr_PC',ChAT);

%% ID plot

blue = [0 0.6 1];
grey = [0.65 0.65 0.65];
figure;
[nm xout] = hist(ID_NT(ID_NT<500),30);   % one outlier of ID=1030.9 excluded
bar(xout,nm,'BarWidth',1,'FaceColor',grey,...
    'EdgeColor',grey)
hold on
[nm2 xout2] = hist(ID_ChAT(ID_ChAT>20&ID_ChAT<Inf),xout);   % one cell clustered based on light-evoked spike shape excluded
bar(xout2,nm2,'BarWidth',1,'FaceColor',blue,...
    'EdgeColor',blue)

edges = 0:0.01:500+0.01;   % bin edges
mx = max(nm);
% mx = 500;
dist_NT = histc(ID_NT(ID_NT<500),edges)';   % histogram
dist_NT = [0 dist_NT(1:end-1)];   % values corresponding to the edges
dist_NT = dist_NT / sum(dist_NT);   % normalize
stairs(edges,cumsum(dist_NT)*mx,'Color',grey/1.5,'LineWidth',2)

dist_ChAT = histc(ID_ChAT(ID_ChAT>20&ID_ChAT<Inf),edges)';   % histogram
dist_ChAT = [0 dist_ChAT(1:end-1)];   % values corresponding to the edges
dist_ChAT = dist_ChAT / sum(dist_ChAT);   % normalize
stairs(edges,cumsum(dist_ChAT)*mx,'Color',blue,'LineWidth',2)

xlabel('Isolation Distance')
ylabel('No. of neurons')

axis([0 500 0 500])
axis square
set(gca,'XTick',[0 250 500],'YTick',[0 250 500])
line([0 500],[mx/2 mx/2],'LineStyle',':','Color','k')
setmyplot_balazs

%% median ID

median([ID_ChAT(ID_ChAT<Inf); ID_NT(ID_NT<Inf)])
se_of_median([ID_ChAT(ID_ChAT<Inf); ID_NT(ID_NT<Inf)])
median(ID_ChAT(ID_ChAT<Inf))
se_of_median(ID_ChAT(ID_ChAT<Inf))

% ID_ChAT2 = ID_ChAT;
% ID_ChAT2(ID_ChAT2==Inf) = 1000;
% median([ID_ChAT2; ID_NT(ID_NT<Inf)])
% se_of_median([ID_ChAT2; ID_NT(ID_NT<Inf)])
% median(ID_ChAT2)
% se_of_median(ID_ChAT2)

%% L-ratio plot

blue = [0 0.6 1];
grey = [0.65 0.65 0.65];
figure;
[nm xout] = hist(Lr_NT,30);
bar(xout,nm,'BarWidth',1,'FaceColor',grey,...
    'EdgeColor',grey)
hold on
[nm2 xout2] = hist(Lr_ChAT(Lr_ChAT<0.15),xout);   % one cell clustered based on light-evoked spike shape excluded
bar(xout2,nm2,'BarWidth',1,'FaceColor',blue,...
    'EdgeColor',blue)

edges = 0:0.0001:0.15+0.01;   % bin edges
mx = max(nm);
% mx = 300;
dist_NT = histc(Lr_NT,edges)';   % histogram
dist_NT = [0 dist_NT(1:end-1)];   % values corresponding to the edges
dist_NT = dist_NT / sum(dist_NT);   % normalize
stairs(edges,cumsum(dist_NT)*mx,'Color',grey/1.5,'LineWidth',2)

dist_ChAT = histc(Lr_ChAT(Lr_ChAT<0.15),edges)';   % histogram
dist_ChAT = [0 dist_ChAT(1:end-1)];   % values corresponding to the edges
dist_ChAT = dist_ChAT / sum(dist_ChAT);   % normalize
stairs(edges,cumsum(dist_ChAT)*mx,'Color',blue,'LineWidth',2)

xlabel('L-ratio')
ylabel('No. of neurons')

axis([0 0.15 0 325])
axis square
set(gca,'XTick',[0 0.15/2 0.15],'YTick',[0 150 300])
line([0 0.15],[mx/2 mx/2],'LineStyle',':','Color','k')
setmyplot_balazs

%% median L-ratio

median([Lr_ChAT; Lr_NT])
se_of_median([Lr_ChAT; Lr_NT])
median(Lr_ChAT)
se_of_median(Lr_ChAT)