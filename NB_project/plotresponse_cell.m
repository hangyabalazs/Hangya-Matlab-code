%%

activation_peak = [allstats.activation_peak];
activation_start = [allstats.activation_start];
activation_end = [allstats.activation_end];
maxvalue = [allstats.maxvalue];
Wpa = [allstats.Wpa];
inhibition_peak = [allstats.inhibition_peak];
inhibition_start = [allstats.inhibition_start];
inhibition_end = [allstats.inhibition_end];
minvalue = [allstats.minvalue];
Wpi = [allstats.Wpi];
baseline = [allstats.baseline];

%% groups

inx_act = Wpa < 0.01;
inx_inh = Wpi < 0.01;
activated = find(inx_act&~inx_inh);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];
inhibited = [inhibited'; inhibited_activated'];

nbareas = {'GP','GP/SI','SI','IC','RT/IC','EP','EA','EAC'};   % aeas considered part of the basal forebrain (BF)
% I = ismember(area1,nbareas);  % if 'primary' area is part of the BF
selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
PV = selectcell(selstr);
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);

[nms inxa PVinx] = intersect(PV,tags);
[nms inxa ChATinx] = intersect(ChAT,tags);
PVactinx = intersect(PVinx,activated);
ChATactinx = intersect(ChATinx,activated);

% tagged = logical(vldty) & (isact==2);
% % inx_act = logical(vldty) & (isact==1) & (baseline>1);
% inx_act = logical(vldty) & (isact==1);
% inx_inh = logical(vldty) & (isinh) & (baseline>1);
% activated = setdiff(find(inx_act&~inx_inh),tagged);
% inhibited = find(inx_inh&~inx_act);
% ai = find(inx_act&inx_inh);
% inx = activation_peak(ai) < inhibition_peak(ai);
% activated_inhibited = sort(ai(inx));
% inhibited_activated = ai(~inx);
% activated = [activated; activated_inhibited];
% inhibited = [inhibited; inhibited_activated];
% tagged = find(tagged);

%% colors

brown = [0.32 0.19 0.19];
purple = [0.48 0.06 0.89];
grey1 = [0.7 0.7 0.7];
grey2 = [0.4 0.4 0.4];
green = [0 0.8 0];
red = [0.8 0 0];
blue = [0 0 0.8];


%%

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',8)
hold on
finh = minvalue(inhibited) ./ baseline(inhibited);
finh(finh==0) = 0.002;
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',grey2,'MarkerFaceColor',grey2,'MarkerSize',8)
% semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',grey1)
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',grey2)
% end
semilogy(activation_peak(ChATactinx),maxvalue(ChATactinx)./baseline(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
semilogy(activation_peak(PVactinx),maxvalue(PVactinx)./baseline(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',[0.32 0.19 0.19])
% end


% semilogy(activation_peak,maxvalue./baseline,'o','MarkerEdgeColor',[0.32 0.19 0.19],'MarkerFaceColor',[0.32 0.19 0.19],'MarkerSize',8)
% for k = 1:length(activation_peak)
%     line([activation_start(k) activation_end(k)],...
%         [maxvalue(k)./baseline(k) maxvalue(k)./baseline(k)],'Color',[0.32 0.19 0.19])
% end
% hold on
% finh = minvalue ./ baseline;
% finh(finh==0) = 0.002;
% semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)

% for k = 1:length(tagged)
%     line([activation_start(tagged(k)) activation_end(tagged(k))],...
%         [maxvalue(tagged(k))./baseline(tagged(k)) maxvalue(tagged(k))./baseline(tagged(k))],'Color',[0 0.7 0])
% end
% line([0 200],[1 1],'Color','k')
% set(gca,'box','off','XTick',0:50:200,'YTick',[0.001 0.01 0.1 1 10 100 1000],'YMinorTick','off',...
%     'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0.01 100])
xlabel('Peak time')
ylabel('Relative change in firing rate')

%%

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',5)
hold on
finh = minvalue(inhibited) ./ baseline(inhibited);
finh(finh==0) = 0.002;
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',grey2,'MarkerFaceColor',grey2,'MarkerSize',5)
% semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',grey1)
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',grey2)
% end
semilogy(activation_peak(ChATactinx),maxvalue(ChATactinx)./baseline(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
semilogy(activation_peak(PVactinx),maxvalue(PVactinx)./baseline(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)
for k = 1:length(ChATactinx)
    line([activation_start(ChATactinx(k)) activation_end(ChATactinx(k))],...
        [maxvalue(ChATactinx(k))./baseline(ChATactinx(k)) maxvalue(ChATactinx(k))./baseline(ChATactinx(k))],'Color',green)
end
for k = 1:length(PVactinx)
    line([activation_start(PVactinx(k)) activation_end(PVactinx(k))],...
        [maxvalue(PVactinx(k))./baseline(PVactinx(k)) maxvalue(PVactinx(k))./baseline(PVactinx(k))],'Color',red)
end

%%

activation_mid = nanmean([activation_start' activation_end'],2)';
inhibition_mid = nanmean([inhibition_start' inhibition_end'],2)';

figure
semilogy(activation_mid(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',5)
hold on
finh = minvalue(inhibited) ./ baseline(inhibited);
finh(finh==0) = 0.002;
semilogy(inhibition_mid(inhibited),finh,'o','MarkerEdgeColor',grey2,'MarkerFaceColor',grey2,'MarkerSize',5)
% semilogy(activation_mid(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',grey1)
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',grey2)
% end
semilogy(activation_mid(pChATactinx),maxvalue(pChATactinx)./baseline(pChATactinx),'o','MarkerEdgeColor',blue,'MarkerFaceColor',blue,'MarkerSize',8)
semilogy(activation_mid(ChATactinx),maxvalue(ChATactinx)./baseline(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
semilogy(activation_mid(PVactinx),maxvalue(PVactinx)./baseline(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)
for k = 1:length(ChATactinx)
    line([activation_start(ChATactinx(k)) activation_end(ChATactinx(k))],...
        [maxvalue(ChATactinx(k))./baseline(ChATactinx(k)) maxvalue(ChATactinx(k))./baseline(ChATactinx(k))],'Color',green)
end
for k = 1:length(PVactinx)
    line([activation_start(PVactinx(k)) activation_end(PVactinx(k))],...
        [maxvalue(PVactinx(k))./baseline(PVactinx(k)) maxvalue(PVactinx(k))./baseline(PVactinx(k))],'Color',red)
end


%%

figure
plot(baseline(activated),activation_peak(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',2)
hold on

plot(baseline(ChATactinx),activation_peak(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot(baseline(PVactinx),activation_peak(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)


%%

figure
plot3(baseline(activated),activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',2)
hold on
grid on
plot3(baseline(ChATactinx),activation_peak(ChATactinx),maxvalue(ChATactinx)./baseline(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot3(baseline(PVactinx),activation_peak(PVactinx),maxvalue(PVactinx)./baseline(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)


%%

figure
plot(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',5)
hold on
plot(inhibition_peak(inhibited),minvalue(inhibited)./baseline(inhibited),'o','MarkerEdgeColor',grey2,'MarkerFaceColor',grey2,'MarkerSize',5)
% plot(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',grey1)
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',grey2)
% end
plot(activation_peak(ChATactinx),maxvalue(ChATactinx)./baseline(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot(activation_peak(PVactinx),maxvalue(PVactinx)./baseline(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)
for k = 1:length(ChATactinx)
    line([activation_start(ChATactinx(k)) activation_end(ChATactinx(k))],...
        [maxvalue(ChATactinx(k))./baseline(ChATactinx(k)) maxvalue(ChATactinx(k))./baseline(ChATactinx(k))],'Color',green)
end
for k = 1:length(PVactinx)
    line([activation_start(PVactinx(k)) activation_end(PVactinx(k))],...
        [maxvalue(PVactinx(k))./baseline(PVactinx(k)) maxvalue(PVactinx(k))./baseline(PVactinx(k))],'Color',red)
end


%%

activation_time = activation_end - activation_start;
inhibition_time = inhibition_end - inhibition_start;

figure
plot(baseline(activated),activation_time(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',2)
hold on

plot(baseline(ChATactinx),activation_time(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot(baseline(PVactinx),activation_time(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)

%%

figure
plot(baseline(activated),activation_peak(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',2)
hold on

plot(baseline(ChATactinx),activation_peak(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot(baseline(PVactinx),activation_peak(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)


%%

figure
plot3(baseline(activated),activation_peak(activated),activation_time(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',2)
hold on
grid on
plot3(baseline(ChATactinx),activation_peak(ChATactinx),activation_time(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot3(baseline(PVactinx),activation_peak(PVactinx),activation_time(PVactinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)

%%

figure
plot3(baseline(activated),activation_time(activated),activation_peak(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',2)
hold on
grid on
plot3(baseline(ChATactinx),activation_time(ChATactinx),activation_peak(ChATactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
plot3(baseline(PVactinx),activation_time(PVactinx),activation_peak(PVctinx),'o','MarkerEdgeColor',red,'MarkerFaceColor',red,'MarkerSize',8)

%%

ainx = find(activation_peak>0.018&activation_peak<0.032);
ainx2 = intersect(ainx,activated);
ainx3 = setdiff(ainx2,ChATinx);
ainx4 = tags(ainx3);
for k = 1:length(ainx4)
    cid = ainx4{k};
    uiopen(['c:\Balazs\_analysis\NB\responsesorter3\' ...
        cid(1:end-2) '_' cid(end) '_PSTH.fig'],1)
end

pChAT = ainx4([30 33 35 36 44 45]);
pChATactinx = ainx3([30 33 35 36 44 45]);
maybe = ainx4([15 22 42]);

%%

% pChAT = 'n029_120211a_3.2' 'n029_120215a_6.1' 'n029_120215a_3.4' 'n029_120216b_6.3'
% 'n029_120222b_4.1' 'n029_120302a_3.3'

% from LDA
pChAT = {'n029_120210a_3.3' ...  % 'n029_120214b_2.7': same as one of the ChAT+ cells; 'n029_120215a_2.2': one of the ChAT+ cells recorded a day later
                'n029_120215a_3.4' 'n029_120221b_6.1' 'n029_120222b_4.1'};   % cellIDs of putative ChAT+ cells
            
% with duplicates (from LDA):
pChAT = {'n029_120202a_3.5' 'n029_120203a_3.1' 'n029_120214b_2.2' 'n029_120210a_3.3'...
    'n029_120214b_2.7' 'n029_120215a_2.2' 'n029_120215a_3.4' 'n029_120221b_6.1'...
    'n029_120222b_4.1'};
pChATactinx = [259 263 427 384 428 461 464 543 559];

%%


hold on
semilogy(activation_peak(pChATactinx),maxvalue(pChATactinx)./baseline(pChATactinx),'o','MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'MarkerSize',8)
for k = 1:length(pChATactinx)
    line([activation_start(pChATactinx(k)) activation_end(pChATactinx(k))],...
        [maxvalue(pChATactinx(k))./baseline(pChATactinx(k)) maxvalue(pChATactinx(k))./baseline(pChATactinx(k))],'Color',[0 0.5 0])
end