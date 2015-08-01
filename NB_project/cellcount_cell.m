%%

selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells

%% cells in behavior

selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells (n = 22)
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative (n = 22)
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells (n = 1316)

%%

for k = 4:19
    tmp = selectcell(['"ratID"==' num2str(k)]);
    nms = cellfun(@(s)s(1:4),tmp,'uniformoutput',false);
    mc = unique(nms);
    lNT = length(intersect(NT,tmp));
    lChAT = length(intersect(ChAT,tmp));
    lpChAT = length(intersect(pChAT,tmp));
    str = [num2str(k) mc num2str(lNT) num2str(lChAT) num2str(lpChAT)];
    disp(str)
end

%% 16 mice

% n013 (ratID==2): 49 NT cells, 0 ChAT, 0 pChAT; ChAT-Cre aguti
% n015 (ratID==4): 15 NT cells, 0 ChAT, 0 pChAT; ChAT-Cre aguti
% n018 (ratID==5): 14 NT cells, 0 ChAT, 1 pChAT; ChAT-Cre aguti
% n020 (ratID==6): 27 NT cells, 0 ChAT, 0 pChAT; ChAT-Cre aguti
% n023 (ratID==8): 159 NT cells, 1 ChAT, 1 pChAT; ChAT-Cre aguti
% n026 (ratID==9): 29 NT cells, 0 ChAT, 0 pChAT; PV-Cre bl6
% n027 (ratID==10): 104 NT cells, 0 ChAT, 0 pChAT; PV-Cre bl6
% n028 (ratID==11): 54 NT cells, 0 ChAT, 2 pChAT; PV-Cre bl6
% n029 (ratID==12): 191 NT cells, 4 ChAT, 8 pChAT; ChAT-Cre bl6
% n037 (ratID==13): 72 NT cells, 0 ChAT, 1 pChAT; ChAT-Cre bl6?
% n038 (ratID==14): 4 NT cells, 0 ChAT, 0 pChAT; ChAT-Cre bl6?
% n039 (ratID==15): 49 NT cells, 0 ChAT, 0 pChAT; ChAT-Cre bl6?
% n040 (ratID==16): 29 NT cells, 0 ChAT, 0 pChAT; PV-Cre bl6
% n043 (ratID==17): 13 NT cells, 0 ChAT, 0 pChAT; ChAT-ChR2 bl6
% n045 (ratID==18): 67 NT cells, 1 ChAT, 1 pChAT; ChAT-ChR2 bl6
% n046 (ratID==19): 210 NT cells, 9 ChAT, 8 pChAT; ChAT-ChR2 bl6

%% training

% n013 (ratID==2): 2-days pre-training
% n015 (ratID==4): 6-days pre-training
% n018 (ratID==5): 2-days pre-training
% n020 (ratID==6): 7-days pre-training
% n023 (ratID==8): 6-days pre-training
% n026 (ratID==9): 3-days pre-training
% n027 (ratID==10): 3-days pre-training
% n028 (ratID==11): 2-days pre-training
% n029 (ratID==12): 1-day pre-training
% n037 (ratID==13): 3-days pre-training
% n038 (ratID==14): 4-days pre-training
% n039 (ratID==15): 5-days pre-training
% n040 (ratID==16): 4-days pre-training
% n043 (ratID==17): 4-days pre-training
% n045 (ratID==18): 4-days pre-training
% n046 (ratID==19): 3-days pre-training
pretraining = [2 6 2 7 6 3 3 2 1 3 4 5 4 4 4 3];

% n013 (ratID==2): above chance by 2nd day of full task
% n015 (ratID==4): above chance by 1st day of full task
% n018 (ratID==5): above chance by 2nd day of full task
% n020 (ratID==6): above chance by 2nd day of full task
% n023 (ratID==8): above chance by 1st day of full task
% n026 (ratID==9): above chance by 1st day of full task
% n027 (ratID==10): above chance by 1st day of full task
% n028 (ratID==11): above chance by 2nd day of full task
% n029 (ratID==12): above chance by 2nd day of full task
% n037 (ratID==13): above chance by 1st day of full task
% n038 (ratID==14): above chance by 1st day of full task
% n039 (ratID==15): above chance by 1st day of full task
% n040 (ratID==16): above chance by 1st day of full task
% n043 (ratID==17): above chance by 1st day of full task
% n045 (ratID==18): above chance by 1st day of full task
% n046 (ratID==19): above chance by 1st day of full task

% n013 (ratID==2): introduce tones from ?th day of full task (no data)
% n015 (ratID==4): introduce tones from 4th day of full task
% n018 (ratID==5): introduce tones from 11th day of full task
% n020 (ratID==6): introduce tones from 5th day of full task
% n023 (ratID==8): introduce tones from 2nd day of full task
% n026 (ratID==9): introduce tones from 4th day of full task
% n027 (ratID==10): introduce tones from 7th day of full task
% n028 (ratID==11): introduce tones from 5th day of full task
% n029 (ratID==12): introduce tones from 5th day of full task
% n037 (ratID==13): introduce tones from 6th day of full task
% n038 (ratID==14): introduce tones from 4th day of full task
% n039 (ratID==15): introduce tones from 5th day of full task
% n040 (ratID==16): introduce tones from 5th day of full task
% n043 (ratID==17): introduce tones from 1st day of full task
% n045 (ratID==18): introduce tones from 1st day of full task
% n046 (ratID==19): introduce tones from 1st day of full task

% n013 (ratID==2): full psychometric from ?th day of full task (animal transferred to other rig, cannot evaluate)
% n015 (ratID==4): full psychometric from 9th day of full task
% n018 (ratID==5): full psychometric from 14th day of full task
% n020 (ratID==6): full psychometric from 9th day of full task
% n023 (ratID==8): full psychometric from 13th day of full task
% n026 (ratID==9): full psychometric from 7th day of full task
% n027 (ratID==10): full psychometric from 10th day of full task
% n028 (ratID==11): full psychometric from 7th day of full task
% n029 (ratID==12): full psychometric from 7th day of full task
% n037 (ratID==13): full psychometric from 10th day of full task
% n038 (ratID==14): full psychometric from 8th day of full task
% n039 (ratID==15): full psychometric from 17th day of full task
% n040 (ratID==16): full psychometric from 8th day of full task
% n043 (ratID==17): full psychometric from 8th day of full task
% n045 (ratID==18): full psychometric from 8th day of full task
% n046 (ratID==19): full psychometric from 10th day of full task
fullpsy = [NaN 9 14 9 13 7 10 7 7 10 8 17 8 8 8 10];

% n013 (ratID==2): max noise 60
% n015 (ratID==4): max noise 60
% n018 (ratID==5): max noise 30
% n020 (ratID==6): max noise 45
% n023 (ratID==8): max noise 20
% n026 (ratID==9): max noise 30
% n027 (ratID==10): max noise 50
% n028 (ratID==11): max noise 50
% n029 (ratID==12): max noise 50
% n037 (ratID==13): max noise 45
% n038 (ratID==14): max noise 45
% n039 (ratID==15): max noise 25
% n040 (ratID==16): max noise 45
% n043 (ratID==17): max noise 45
% n045 (ratID==18): max noise 45
% n046 (ratID==19): max noise 45
maxnoise = [60 60 30 45 20 30 50 50 50 45 45 25 45 45 45 45];

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

edges = 0:0.0001:500+0.01;   % bin edges
mx = max(nm);
dist_NT = histc(ID_NT(ID_NT<500),edges)';   % histogram
dist_NT = [0 dist_NT(1:end-1)];   % values corresponding to the edges
dist_NT = dist_NT / sum(dist_NT);   % normalize
stairs(edges,cumsum(dist_NT)*mx,'Color',grey/1.5,'LineWidth',2)

dist_ChAT = histc(ID_ChAT(ID_ChAT>20&ID_ChAT<Inf),edges)';   % histogram
dist_ChAT = [0 dist_ChAT(1:end-1)];   % values corresponding to the edges
dist_ChAT = dist_ChAT / sum(dist_ChAT);   % normalize
stairs(edges,cumsum(dist_ChAT)*mx,'Color',blue,'LineWidth',2)

%% cells in behavior - HDB

choosecb('HDB')
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified (n = 12)
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
NT = selectcell(['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % NT (n = 209)