%% Reaction time

inx1 = (TE.ResponseType==1);   % hit
inx2 = (TE.ResponseType==2);   % false alarm
inx3 = (TE.ResponseType==3);   % correct reject
inx4 = (TE.ResponseType==4);   % miss

edges = min(TE.ReactionTime):0.1:max(TE.ReactionTime);
cnts = (edges(1:end-1) + edges(2:end)) / 2;

rt1 = histc(TE.ReactionTime(inx1),edges);
rt2 = histc(TE.ReactionTime(inx2),edges);
rt3 = histc(TE.ReactionTime(inx3),edges);
rt4 = histc(TE.ReactionTime(inx4),edges);
rt1 = rt1(1:end-1);
rt2 = rt2(1:end-1);
rt3 = rt3(1:end-1);
rt4 = rt4(1:end-1);

figure
plot(cnts,rt1,'g')
hold on
plot(cnts,rt2,'r')
% plot(cnts,rt3,'Color',[0 0.25 0])
% plot(cnts,rt4,'Color',[0.25 0 0])

%% Reaction time for different intensities

soundints = sort(unique(TE.StimulusDuration),'descend');
lensi = length(soundints);
si = cell(1,lensi);
rt = cell(1,lensi);
clr = zeros(lensi,3);
edges = min(TE.ReactionTime):0.5:max(TE.ReactionTime);
cnts = (edges(1:end-1) + edges(2:end)) / 2;
for k = 1:lensi
    si{k} = (TE.StimulusDuration==soundints(k)&TE.ResponseType==2);
%     si{k} = (TE.ResponseType==2);
    clr(k,:) = [1-(k-1)/lensi (k-1)/lensi 0];
    rt{k} = histc(TE.ReactionTime(si{k}),edges);
    rt{k} = rt{k}(1:end-1);
end

figure
plot(cnts,rt{1},'Color',clr(1,:))
hold on
for k = 2:lensi
    plot(cnts,rt{k},'Color',clr(k,:))
end


    