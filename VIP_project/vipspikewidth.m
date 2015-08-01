%% load

load('C:\Balazs\_analysis\VIP\vipisinfluenced11\validity.mat')
load('C:\Balazs\_analysis\VIP\vipisinfluenced11\p_val.mat')
load('C:\Balazs\_analysis\VIP\vipisinfluenced11\spike_width.mat')

%% without FR criterion

inx_act = logical(vldty') & (p_act<0.015);
inx_inh = logical(vldty') & (p_inh<0.015);

%% with FR criterion

inx_act = logical(vldty') & (p_act<0.015) & (baseline>2);
inx_act(105) = 1;   % tagged
inx_inh = logical(vldty') & (p_inh<0.015) & (baseline>2);

%% groups

tagged = find(inx_act&activation_start<3);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = ai(inx);
inhibited_activated = ai(~inx);
activated = [activated activated_inhibited];
inhibited = [inhibited inhibited_activated];

%% plots

figure
plot(spike_width(activated),activation_peak(activated),'ro','MarkerfaceColor','r')
hold on
plot(spike_width(inhibited),inhibition_peak(inhibited),'bo','MarkerfaceColor','b')

figure
plot(baseline(activated),activation_peak(activated),'ro','MarkerfaceColor','r')
hold on
plot(baseline(inhibited),inhibition_peak(inhibited),'bo','MarkerfaceColor','b')

%%

figure;
plot(spike_width,baseline,'k.')
hold on
plot(spike_width(activated),baseline(activated),'ro','MarkerfaceColor','r')
plot(spike_width(inhibited),baseline(inhibited),'bo','MarkerfaceColor','b')
line([100 600],[2 2],'Color','g')
% plot(spike_width(both),baseline(both),'go','MarkerfaceColor','g')