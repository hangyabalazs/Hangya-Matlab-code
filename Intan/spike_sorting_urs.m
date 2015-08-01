%data = read_intan_data;
% sampling rate 25kHz
[b,a] = butter(3,[600 6000]./25000);
filt_data = filter(b,a,data(:,1));
plot(filt_data);
%%
[peak_amplitude, peak_index] = max(WaveForms,[],3);
[min_amplitude, min_index] = min(WaveForms,[],3);
delta_peak = min_index - peak_index;
area = sum(abs(WaveForms),3);
[COEFF,SCORE] = princomp(WaveForms(:,1,:));
%%
tmp = ypos;
for i = 1:size(tmp,1)
    if tmp(i,1) < 0
        tmp(i,:) = tmp(i-1,:);
    end
end
ypos_clean2 = tmp;
%%
for i =1:size(TS,1)
    idx(i)=find(TimeStamps==TS(i));
end
%%
for i = 1:length(xytms)
    [c,index(i)]=min(abs(TTLs-xytms(i)));
end
%%
for i = 1:length(TS)
    [c,spike_open(i)]=min(abs(xytms_open-TS(i)));
end
%%
heat_map = zeros(max(xy_open(spike_open,1)),max(xy_open(spike_open,2)));
for i = 1:length(spike_open)
    heat_map(xy_open(spike_open(i),1),xy_open(spike_open(i),2)) = heat_map(xy_open(spike_open(i),1),xy_open(spike_open(i),2))+1;
end
%%