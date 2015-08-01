function CCG_effect_size
%CCG_EFFECT_SIZE   Relative effect size for cross-correlograms.
%   CCG_EFFECT_SIZE quantifies negative/positive peak value (min/max) in
%   the -4 to +4 ms window relative to the average count outside this
%   window (ratio).
%
%   See also SOM_CCG_CONF_FILTER and CCG_EFFECT_SIZE2.

% Effect size for PV-PV
pv_exc = [1.51 1.18 1.19] - 1;   % excitation
pv_inh = [0.72 0.84 0.87 0.79] - 1;   % inhibition

% Effect size for PV-NT
pvnt_exc = [2.08 7.11 1.38 1.26 1.25 1.11 1.32 1.16 1.15 2.63 1.83 1.28 3.19 ...
    1.13 1.39 1.30 1.11 1.30 1.13 1.17 1.17 2.37 1.15 1.20 1.09 1.14 1.62 ...
    1.31 1.30 1.62 1.56 1.49 2.38] - 1;   % excitation
pvnt_inh = [0.73 0.53 0.62 0.83 0.85 0.74 0.69 0.81 0.82 0.67 0.03 0.44 0.32 ...
    0.05 0.77 0.02 0.44 0.88] - 1;   % inhibition

% Effect size for SOM-NT
somnt_exc = [2.61 3.41 1.15 1.32 3.22 1.72 3.49 3.23 6.59 1.17 1.38 3.59 ...
    6.09 2.19 2.07] - 1;   % excitation
somnt_inh = [0.42 0.46] - 1;   % inhibition

% Plot bar graph with mean and SE
figure
hold on
bar(1,mean(pv_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, excitation
bar(1,mean(pv_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, inhibition
E = errorbar([1 1],[mean(pv_exc) mean(pv_inh)],...
    [std(pv_exc)/sqrt(length(pv_exc)) std(pv_inh)/sqrt(length(pv_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(2,mean(pvnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, excitation
bar(2,mean(pvnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, inhibition
E = errorbar([2 2],[mean(pvnt_exc) mean(pvnt_inh)],...
    [std(pvnt_exc)/sqrt(length(pvnt_exc)) std(pvnt_inh)/sqrt(length(pvnt_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(3,mean(somnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, excitation
bar(3,mean(somnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, inhibition
E = errorbar([3 3],[mean(somnt_exc) mean(somnt_inh)],...
    [std(somnt_exc)/sqrt(length(somnt_exc)) std(somnt_inh)/sqrt(length(somnt_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

% keyboard

% Plot bar graph with median and IQR
figure
hold on
bar(1,median(pv_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, excitation
bar(1,median(pv_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, inhibition
E = errorbar([1 1],[median(pv_exc) median(pv_inh)],...
    [prctile(pv_exc,25) prctile(pv_inh,25)],...
    [prctile(pv_exc,75) prctile(pv_inh,75)],...
    '+','Color',[0 0 0],'LineWidth',2);   % IQR
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(2,median(pvnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, excitation
bar(2,median(pvnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, inhibition
E = errorbar([2 2],[median(pvnt_exc) median(pvnt_inh)],...
    [prctile(pvnt_exc,25) prctile(pvnt_inh,25)],...
    [prctile(pvnt_exc,75) prctile(pvnt_inh,75)],...
    '+','Color',[0 0 0],'LineWidth',2);   % IQR
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(3,median(somnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, excitation
bar(3,median(somnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, inhibition
E = errorbar([3 3],[median(somnt_exc) median(somnt_inh)],...
    [prctile(somnt_exc,25) prctile(somnt_inh,25)],...
    [prctile(somnt_exc,75) prctile(somnt_inh,75)],...
    '+','Color',[0 0 0],'LineWidth',2);   % IQR
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

keyboard