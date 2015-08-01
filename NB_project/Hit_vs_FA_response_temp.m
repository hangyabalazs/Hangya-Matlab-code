function Hit_vs_FA_response
%HIT_VS_FA_RESPONSE   Comparisan of neural response for positive and negative reinforcement.
%   HIT_VS_FA_RESPONSE compares spike reliability, latency and jitter for
%   hit and false alarm trials. The analysis is restricted to those
%   (identified and putative) cholinergic cells that are significantly
%   activated by both negative and positive reinforcement (p<0.01,
%   one-sided Mann-Whitney U-test).
%
%   See also NBACTIVATION and RELIABILITY_LATENCY_JITTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   23-July-2013

%   Edit log: BH 7/23/13

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT pChAT];   % identified and putative

% Response properties
FA_latency = getvalue('FA_Latency',allChAT);   % latency for false alarms
Hit_latency = getvalue('Hit_Latency',allChAT);   % latency for hits
FA_jitter = getvalue('FA_Jitter',allChAT);   % jitter for hits
Hit_jitter = getvalue('Hit_Jitter',allChAT);   % jitter for false alarms
FA_reliability = getvalue('FA_Reliability',allChAT);   % reliability for false alarms
Hit_reliability = getvalue('Hit_Reliability',allChAT);   % reliability for hits

% PSTH statistics
FA_psth_stats = getvalue('FA_psth_stats',allChAT);
Hit_psth_stats = getvalue('Hit_psth_stats',allChAT);
FA_psth_stats = nancell2struct(FA_psth_stats);
Hit_psth_stats = nancell2struct(Hit_psth_stats);
FA_act = [FA_psth_stats.Wpa] < 0.01;  % cells activated after false alarms
Hit_act = [Hit_psth_stats.Wpa] < 0.01;   % cells activated after hits

% Scatter plots
inx = FA_act & Hit_act;  % cells activated in both hit and false alarm trials
figure
plot(Hit_latency(inx),FA_latency(inx),'o',...
    'MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]')  % Hit vs. FA latency
xlabel('latency, hits')
ylabel('latency, false alarms')
axis([0.01 0.08 0.01 0.032])
line([0 0.032],[0 0.032],'Color','k','LineWidth',3)

figure
plot(Hit_jitter(inx),FA_jitter(inx),'o',...
    'MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]')  % Hit vs. FA jitter
xlabel('jitter, hits')
ylabel('jitter, false alarms')
axis([0 0.01 0 0.01])
line([0 0.01],[0 0.01],'Color','k','LineWidth',3)

figure
plot(Hit_reliability(inx),FA_reliability(inx),'o',...
    'MarkerSize',12,'MarkerFaceColor',[0 0.8 0],'MarkerEdgeColor',[0 0.8 0]')  % Hit vs. FA reliability
xlabel('reliability, hits')
ylabel('reliability, false alarms')
axis([0.2 1 0.2 1])
line([0.2 1],[0.2 1],'Color','k','LineWidth',3)

% Statistics
boxstat(Hit_latency(inx),FA_latency(inx),'Hit','FalseAlarm',[],'paired')  % compare latencies
title('Latency')

boxstat(Hit_jitter(inx),FA_jitter(inx),'Hit','FalseAlarm',[],'paired')  % compare jitters
title('Jitter')

boxstat(Hit_reliability(inx),FA_reliability(inx),'Hit','FalseAlarm',[],'paired')  % compare reliabilities
title('Reliability')

% Regression
y = Hit_latency(inx)';  % latency
x = FA_latency(inx)';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance

y = Hit_latency(inx&Hit_latency'<0.05)';  % latency, without the 3 outliers
x = FA_latency(inx&Hit_latency'<0.05)';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance

y = Hit_jitter(inx)';  % jitter
x = FA_jitter(inx)';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance

y = Hit_reliability(inx)';  % reliability
x = FA_reliability(inx)';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance

keyboard