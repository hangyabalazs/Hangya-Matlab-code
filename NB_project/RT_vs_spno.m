function RT_vs_spno
%RT_VS_SPNO   Association between reaction time or stimulus intensity and number of evoked spikes
%   RT_VS_SPNO calculates regresseion between reaction time or stimulus
%   intensity and number of spikes in bursts evoked by water or air puff in
%   cholinergic neurons. (identified or putative).
%
%   See also FA_VS_HIT_RESPONSE.

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT pChAT];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Number of spikes fired in response to reinforcement
FA_spno = getvalue('FA_SpikeNumberDistribution',allChAT);   % spike numbers for false alarms
Hit_spno = getvalue('Hit_SpikeNumberDistribution',allChAT);   % spike numbers for false alarms

% PSTH statistics
FA_psth_stats = getvalue('FA_psth_stats',allChAT);   % PSTH statistics, false alarms
Hit_psth_stats = getvalue('Hit_psth_stats',allChAT);   % PSTH statistics, hits
FA_psth_stats = nancell2struct(FA_psth_stats);
Hit_psth_stats = nancell2struct(Hit_psth_stats);
FA_act = [FA_psth_stats.Wpa] < 0.01;  % cells activated after false alarms
Hit_act = [Hit_psth_stats.Wpa] < 0.01;   % cells activated after hits
FA_act_inx = find(FA_act);
Hit_act_inx = find(Hit_act);
ChAT_FA_act = allChAT(FA_act);   % cholinergic cells activated by false alarms
ChAT_Hit_act = allChAT(Hit_act);   % cholinergic cells activated by hits
NumFA_act = length(ChAT_FA_act);   % number of cells activated by false alarms
NumHit_act = length(ChAT_Hit_act);   % number of cells activated by hits

% Association between RT/stimulus intensity and cholinergic response, False Alarms
ngallp = nan(NumFA_act,1);
for iC = 1:NumFA_act  % loop through FA-activated cholinergic cells
    inx = FA_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    TE = loadcb(cellid,'TrialEvents');  % trial events
    ngRT = TE.NoGoRT;   % no-go reaction time
    ngRT = ngRT(~isnan(ngRT));
    ngSI = TE.SoundIntensity;   % stim. intensity
    ngSI = ngSI(~isnan(TE.FalseAlarm));
    ngspno = FA_spno{inx};   % number of spikes fired in response to reinforcement
%     p = plotassoc(ngSI,ngspno);  % plot association between stim. int. and spike number and return p-value
%     ngallp(inx) = p;  % p-values for regression
    p = plotassoc(ngspno,ngRT);  % plot association between RT and spike number and return p-value
    ngallp(inx) = p;  % p-values for regression
end

% Association between RT/stimulus intensity and cholinergic response, Hits
goallp = nan(NumHit_act,1);
for iC = 1:NumHit_act  % loop through Hit-activated cholinergic cells
    inx = Hit_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    TE = loadcb(cellid,'TrialEvents');  % trial events
    goRT = TE.GoRT;   % go reaction time
    goRT = goRT(~isnan(goRT));
    goSI = TE.SoundIntensity;   % stim. intensity
    goSI = goSI(~isnan(TE.Hit));
    gospno = Hit_spno{inx};   % number of spikes fired in response to reinforcement
%     p = plotassoc(goSI,gospno);  % plot association between stim. int. and spike number and return p-value
%     goallp(inx) = p;  % p-values for regression
    p = plotassoc(gospno,goRT);  % plot association between RT and spike number and return p-value
    goallp(inx) = p;  % p-values for regression
end

% -------------------------------------------------------------------------
function p = plotassoc(X,Y)

% inx1 = Y==0;   % no response
% inx2 = Y>=1;   % response
%     if sum(inx1) * sum(inx2) > 0   % neither of the index sets are empty
%         boxstat(ngRT(inx1),ngRT(inx2),'no response','response');  % stat., false alarms
%         title([regexprep(cellid,'_',' ') ' FA'])
%     end

% Number of eveoked spikes vs. stimulus intensity
unq = sort(unique(X),'ascend');   % different stim. intensities
NumSI = length(unq);   % number of different stim. intensities
[xk yk] = deal([]);
[ngmn ngse] = deal(nan(1,NumSI));
for k = 1:NumSI  % loop through different stim. intensities
    sic = unq(k);   % stim. intensity
    ngss = Y(X==sic);   % number of evoked spikes for the selected stim. intensity
    nm = sum(X==sic);
    xk = [xk ngss];   % #evoked spikes
    yk = [yk ones(1,nm)*k];   % grouping variable
    ngmn(k) = mean(ngss);   % mean spike number
    ngse(k) = std(ngss) / sqrt(nm);   % SE of spike number
end

% Plot
figure
boxplot(xk,yk)
figure
bar(unq,ngmn,0.5,'FaceColor','none')
hold on
errorbar(unq,ngmn,ngse,'k+')

% Regression
y = Y;   % number of evoked spikes
x = X;   % stimulus intensity/RT
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
disp(p)

% if p < 0.05
%     plot(unq,ngmn,'k')
% else
%     plot(unq,ngmn,'Color',[0.7 0.7 0.7])
% end

% keyboard