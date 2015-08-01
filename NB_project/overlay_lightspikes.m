%%

fn = 'f:\VIP_gonogo_cellbase2\ml01\130506a\'
[EventTimeStamps, EventIDs, Nttls, Extras, EventStrings NlxHeader] = ...
    Nlx2MatEV([fn 'Events.nev'],[1 1 1 1 1],1,1,1);

%%

inx = find(Nttls==16384);
light_ts = EventTimeStamps(inx) / 10^2;

%%

global MClust_FeatureTimestamps
f_ts = MClust_FeatureTimestamps;

lag1 = 0;
lag2 = 80;   % 4 ms
show_ts = [];
for k = 1:length(light_ts)
    show_ts = [show_ts; find(f_ts>light_ts(k)+lag1&f_ts<light_ts(k)+lag2)];
end

%%

global MClust_CurrentFeatureData
figure
h = plot(MClust_CurrentFeatureData(show_ts,1), MClust_CurrentFeatureData(show_ts,2), 'b+');

%%

hold on
h = plot(MClust_CurrentFeatureData(show_ts,1), MClust_CurrentFeatureData(show_ts,2), 'y+');