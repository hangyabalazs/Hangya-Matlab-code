%%

[data1, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH1.continuous');
[data2, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH2.continuous');
[data3, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH3.continuous');
[data4, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH4.continuous');
data = [data1 data2 data3 data4];

%%

[data1, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH5.continuous');
[data2, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH6.continuous');
[data3, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH7.continuous');
[data4, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\gezo_2015-06-07_18-12-18\100_CH8.continuous');
data = [data1 data2 data3 data4];

%%

[data1, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\frank_2015-06-07_17-51-01\100_CH10.continuous');
[data2, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\frank_2015-06-07_17-51-01\100_CH11.continuous');
[data3, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\frank_2015-06-07_17-51-01\100_CH12.continuous');
[data4, timestamps, info] = load_open_ephys_data('d:\TENSS2015\data\frank_2015-06-07_17-51-01\100_CH13.continuous');
data = [data1 data2 data3 data4];

%%

avw = squeeze(mean(AllWaveForms{1},1));
figure;plot(avw')

%%

mx = squeeze(max(AllWaveForms{1},[],3));
mn = squeeze(min(AllWaveForms{1},[],3));
m = mx + mn;
figure;plot(m(:,1),mx(:,2),'.')