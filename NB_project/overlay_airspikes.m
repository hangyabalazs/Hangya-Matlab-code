%%

load('G:\junk_cellbase\n045\121217b\AirpuffSpikeTimes.mat')
light_ts = AirpuffSpikeTimes / 10^2 * 10^6;

%%

global MClust_FeatureTimestamps
f_ts = MClust_FeatureTimestamps;

lag1 = -2;
lag2 = 2;   % 0.20 ms
show_ts = [];
for k = 1:length(light_ts)
    show_ts = [show_ts; find(f_ts>light_ts(k)+lag1&f_ts<light_ts(k)+lag2)];
end
global AirpuffSpikes
AirpuffSpikes = show_ts;

%%

global MClust_CurrentFeatureData
figure
h = plot(MClust_CurrentFeatureData(show_ts,1), MClust_CurrentFeatureData(show_ts,2), 'b+');

%%

hold on
h = plot(MClust_CurrentFeatureData(show_ts,1), MClust_CurrentFeatureData(show_ts,2), 'r+');