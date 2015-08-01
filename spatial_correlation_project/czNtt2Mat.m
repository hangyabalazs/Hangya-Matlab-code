%% helping codes for spike sorting (interneuron-place cell pair project)

dr = 'C:\MATLAB_R2007a\work\Nlx_converter';
addpath(genpath(dr))

%%

[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints NlxHeader] = ...
    Nlx2MatSpike('X:\In_Vivo\balazs\_analysis\Czurko2\INTERNEURONS\_osc_int_EEG\acin08s007\2002-9-30_16-37-7\Sc2.Ntt',...
    [1 1 1 1 1],1,1,1);

%%

p1 = 2;
p2 = 7;


inx = find(CellNumbers==2);
length(inx)
figure;plot(Params(p1,:),Params(p2,:),'k.','MarkerSize',4)
hold on
plot(Params(p1,inx),Params(p2,inx),'r.','MarkerSize',4)

inx = find(CellNumbers==4);
length(inx)
plot(Params(p1,inx),Params(p2,inx),'g.','MarkerSize',4)

inx = find(CellNumbers==6);
length(inx)
plot(Params(p1,inx),Params(p2,inx),'b.','MarkerSize',4)
% plot(Params(1,inx),Params(4,inx),'r.','MarkerSize',4)

%%

inx = find(CellNumbers==6);
length(inx)
figure
subplot(1,4,1)
plot(squeeze(DataPoints(:,1,inx)),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,1,inx),3)),'r')
subplot(1,4,2)
plot(squeeze(DataPoints(:,2,inx)),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,2,inx),3)),'r')
subplot(1,4,3)
plot(squeeze(DataPoints(:,3,inx)),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,3,inx),3)),'r')
subplot(1,4,4)
plot(squeeze(DataPoints(:,4,inx)),'Color',[0.8 0.8 0.8])
hold on
plot(squeeze(mean(DataPoints(:,4,inx),3)),'r')


%% PCA

DP1 = squeeze(DataPoints(:,1,:))';
[coeff_DP1, scores_DP1] = princomp(DP1);

DP2 = squeeze(DataPoints(:,2,:))';
[coeff_DP2, scores_DP2] = princomp(DP2);

DP3 = squeeze(DataPoints(:,3,:))';
[coeff_DP3, scores_DP3] = princomp(DP3);

DP4 = squeeze(DataPoints(:,4,:))';
[coeff_DP4, scores_DP4] = princomp(DP4);

%%

close all
S1 = scores_DP1(:,1);
S2 = scores_DP4(:,2);
% S1 = Params(4,:);
% S2 = Params(7,:);

figure
plot(S1,S2,'k.','MarkerSize',4)
hold on

inx = find(CellNumbers==2);
plot(S1(inx),S2(inx),'r.','MarkerSize',4)
inx = find(CellNumbers==4);
plot(S1(inx),S2(inx),'g.','MarkerSize',4)
inx = find(CellNumbers==6);
plot(S1(inx),S2(inx),'b.','MarkerSize',4)

%%

biplot(coeff(:,1:2),'scores',scores(:,1:2))

%% energy

l1 = squeeze(DataPoints(:,1,:));
E1 = sum(l1.^2,1) / size(l1,2);

l2 = squeeze(DataPoints(:,2,:));
E2 = sum(l2.^2,1) / size(l2,2);

l3 = squeeze(DataPoints(:,3,:));
E3 = sum(l3.^2,1) / size(l3,2);

l4 = squeeze(DataPoints(:,4,:));
E4 = sum(l4.^2,1) / size(l4,2);

%% area

l1 = squeeze(DataPoints(:,1,:));
A1 = sum(abs(l1),1);

l2 = squeeze(DataPoints(:,2,:));
A2 = sum(abs(l2),1);

l3 = squeeze(DataPoints(:,3,:));
A3 = sum(abs(l3),1);

l4 = squeeze(DataPoints(:,4,:));
A4 = sum(abs(l4),1);

%% PCA of energy-normalized waveform

DPE1 = squeeze(DataPoints(:,1,:))' ./ repmat(E1',1,32);
[coeff_DPE1, scores_DPE1] = princomp(DPE1);

DPE2 = squeeze(DataPoints(:,2,:))' ./ repmat(E2',1,32);
[coeff_DPE2, scores_DPE2] = princomp(DPE2);

DPE3 = squeeze(DataPoints(:,3,:))' ./ repmat(E3',1,32);
[coeff_DPE3, scores_DPE3] = princomp(DPE3);

DPE4 = squeeze(DataPoints(:,4,:))' ./ repmat(E4',1,32);
[coeff_DPE4, scores_DPE4] = princomp(DPE4);

%% PCA of area-normalized waveform

DPA1 = squeeze(DataPoints(:,1,:))' ./ repmat(A1',1,32);
[coeff_DPA1, scores_DPA1] = princomp(DPA1);

DPA2 = squeeze(DataPoints(:,2,:))' ./ repmat(A2',1,32);
[coeff_DPA2, scores_DPA2] = princomp(DPA2);

DPA3 = squeeze(DataPoints(:,3,:))' ./ repmat(A3',1,32);
[coeff_DPA3, scores_DPA3] = princomp(DPA3);

DPA4 = squeeze(DataPoints(:,4,:))' ./ repmat(A4',1,32);
[coeff_DPA4, scores_DPA4] = princomp(DPA4);

%% slice

n = size(DataPoints,3);
z = 1;
slope = zeros(1,n);
for k = 1:n
    wn = squeeze(DataPoints(:,z,k))';
    mxinx = find(wn==max(wn));
    mxinx = mxinx(end);
    lm = 12;
    linlen = 5;
    lsrate = 10000;
    time = linspace(0,linlen/lsrate,linlen+1);
    r = zeros(1,lm);
    gr = zeros(1,lm);
    icp = zeros(1,lm);
    for m = 1:lm
        [gr0,icp0,err0] = linefit(time,wn(mxinx+m-1:mxinx+m+linlen-1));
        r(m) = err0;
        gr(m) = gr0;
        icp(m) = icp0;
    end
    lsr = min(r);
    grinx = find(r==lsr);
    slope(k) = gr(grinx) / 1000;   % slope is given in mV/ms
    
    figure
    plot(wn)
    hold on
    t = [mxinx+grinx:mxinx+grinx+linlen];
    y = time .* slope(k) + intercept(k);
    plot(t,y,'r')
    r
end

%% plot sort

figure

inx = find(strcmp(handles.xfeature,handles.feature_field_names(:,1)));
xfield = handles.feature_field_names{inx,2};
eval(['S1 = handles.' xfield ';']);
inx = find(strcmp(handles.yfeature,handles.feature_field_names(:,1)));
yfield = handles.feature_field_names{inx,2};
eval(['S2 = handles.' yfield ';']);

% Plot
plot(S1,S2,'k.','MarkerSize',4)
hold on

for k = 1:length(handles.cell_codes)
    inx = find(handles.CellNumbers==handles.cell_codes(k));
    plot(S1(inx),S2(inx),handles.plot_codes{k},'MarkerSize',4)
end
xlabel(handles.xfeature)
ylabel(handles.yfeature)
fn = ['sort_' handles.xfeature '_' handles.yfeature];
saveas(gcf,fn)

%% plot sort black

xlabel('')
ylabel('')
ach = allchild(gca);
blk = findobj(ach,'Color','k');
set(blk,'Color',[0.31 0.31 0.31])
set(gca,'Color','k')
pos = [183 107 872 733];
set(gca,'YTick',[]);set(gca,'XTick',[]);global pos;set(gcf,'Position',pos);
fn = ['sort_' handles.xfeature '_' handles.yfeature '_black'];
saveas(gcf,fn)


%% plot waveforms

pos = [440   104   291   687];
alim = 2000;
figure
set(gcf,'Position',pos)
DataPoints = handles.DataPoints;
plot(squeeze(mean(DataPoints(:,1,inx),3)),'k','LineWidth',3)
hold on
plot(squeeze(mean(DataPoints(:,1,inx),3))+squeeze(std(DataPoints(:,1,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,1,inx),3))-squeeze(std(DataPoints(:,1,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,1,inx),3)),'k','LineWidth',3)
xlim([1 32])
ylim([-1*alim alim])
axis off
fn = ['waveform_cell' num2str(cn) '_t1.fig'];
saveas(gcf,fn)

figure
set(gcf,'Position',pos)
DataPoints = handles.DataPoints;
plot(squeeze(mean(DataPoints(:,2,inx),3)),'k','LineWidth',3)
hold on
plot(squeeze(mean(DataPoints(:,2,inx),3))+squeeze(std(DataPoints(:,2,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,2,inx),3))-squeeze(std(DataPoints(:,2,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,2,inx),3)),'k','LineWidth',3)
xlim([1 32])
ylim([-1*alim alim])
axis off
fn = ['waveform_cell' num2str(cn) '_t2.fig'];
saveas(gcf,fn)

figure
set(gcf,'Position',pos)
DataPoints = handles.DataPoints;
plot(squeeze(mean(DataPoints(:,3,inx),3)),'k','LineWidth',3)
hold on
plot(squeeze(mean(DataPoints(:,3,inx),3)),'k','LineWidth',3)
plot(squeeze(mean(DataPoints(:,3,inx),3))+squeeze(std(DataPoints(:,3,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,3,inx),3))-squeeze(std(DataPoints(:,3,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
xlim([1 32])
ylim([-1*alim alim])
axis off
fn = ['waveform_cell' num2str(cn) '_t3.fig'];
saveas(gcf,fn)

figure
set(gcf,'Position',pos)
DataPoints = handles.DataPoints;
plot(squeeze(mean(DataPoints(:,4,inx),3)),'k','LineWidth',3)
hold on
plot(squeeze(mean(DataPoints(:,4,inx),3))+squeeze(std(DataPoints(:,4,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,4,inx),3))-squeeze(std(DataPoints(:,4,inx),[],3)),'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(squeeze(mean(DataPoints(:,4,inx),3)),'k','LineWidth',3)
xlim([1 32])
ylim([-1*alim alim])
axis off
fn = ['waveform_cell' num2str(cn) '_t4.fig'];
saveas(gcf,fn)

%% plot auto

xd=get(gco,'XData');
yd=get(gco,'YData');
yd=yd/sum(yd);
set(gco,'YData',yd)

%%

load('X:\In_Vivo\balazs\_analysis\Czurko\discriminated2\new\acin11s007_spk\nr_03_2.mat')
data = data';
[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints NlxHeader] = ...
    Nlx2MatSpike('X:\In_Vivo\balazs\_analysis\Czurko\Ntt\acin11s007_sc4New.Ntt',...
    [1 1 1 1 1],1,1,1);

inx = find(CellNumbers==2);
TimeStamps2 = TimeStamps - TimeStamps(1);
data2 = TimeStamps2(inx);
min(data-data2/1000000)
max(data-data2/1000000)
corrv = mean(data-data2/1000000);
