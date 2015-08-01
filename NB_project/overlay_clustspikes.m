%%

load('F:\junk_cellbase\n045\121217x\TT4_selected_alldiff.mat')

%%

global MClust_FeatureTimestamps
f_ts = MClust_FeatureTimestamps;

[jn1 jn2 show_ts] = intersect(TS,f_ts);

%%

global MClust_CurrentFeatureData
figure
h = plot(MClust_CurrentFeatureData(show_ts,1), MClust_CurrentFeatureData(show_ts,2), 'b+');

%%

hold on
h = plot(MClust_CurrentFeatureData(show_ts,1), MClust_CurrentFeatureData(show_ts,2), 'y+');