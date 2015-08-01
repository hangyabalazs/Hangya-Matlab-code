function ITI_vs_spno
%ITI_VS_SPNO   Association between foreperiod and number of evoked spikes
%   ITI_VS_SPNO calculates regresseion between foreperiod and number of 
%   spikes in bursts evoked by water or air puff in cholinergic neurons
%   (identified or putative).
%
%   See also RT_VS_SPNO.

% Directories
global DATAPATH
resdir = [DATAPATH 'NB\ITI_vs_spno\'];
issave = true;

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

% Progress indicator
wb = nicewaitbar(0,'Running ''ITI vs spno'' for false alarms...','Name','Please wait...');  % progress indicator
global WB
WB(end+1) = wb;

% Association between foreperiod and cholinergic response, False Alarms
ngallp = nan(NumFA_act,1);
for iC = 1:NumFA_act  % loop through FA-activated cholinergic cells
    inx = FA_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    
    % Foreperiod, reward size, number of evoked spikes
    TE = loadcb(cellid,'TrialEvents');  % trial events
    ngITI = TE.ITIDistribution;   % foreperiod
    ngITI = ngITI(~isnan(TE.FalseAlarm));
    ngTITI = TE.TotalITI;   % realized foreperiod (with restarts)
    ngTITI = ngTITI(~isnan(TE.FalseAlarm));
    ngRew = TE.LeftWaterVolume;   % reward size
    ngRew = ngRew(~isnan(TE.FalseAlarm));
    ngspno = FA_spno{inx};   % number of spikes fired in response to reinforcement
    
    % Plot association between foreperiod and number of evoked spikes and return p-value
    p = plotassoc(ngspno,ngITI);  % plot association between foreperiod and spike number and return p-value
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('Foreperiod')
    xlabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_ITI_VS_NGSPNO'];
        saveas(gcf,fnm)   % save bar graph
    end
    
    % Plot association between full foreperiod and number of evoked spikes and return p-value
    p = plotassoc(ngspno,ngTITI);  % plot association between full foreperiod and spike number and return p-value
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('Full foreperiod')
    xlabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_TITI_VS_NGSPNO'];
        saveas(gcf,fnm)   % save bar graph
    end
    
    % Plot association between foreperiod and number of evoked spikes and return p-value
    p = plotassoc(ngRew,ngspno);  % plot association between reward size and spike number and return p-value
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    xlabel('Reward size')
    ylabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_REW_VS_NGSPNO2'];
        saveas(gcf,fnm)   % save bar graph
    end
    
    ngallp(inx) = p;  % p-values for regression
    close all
    waitbar(iC/NumFA_act)
end
close(wb)

% Progress indicator
wb = nicewaitbar(0,'Running ''ITI vs spno'' for hits...','Name','Please wait...');  % progress indicator
WB(end+1) = wb;

% Association between foreperiod and cholinergic response, Hits
goallp = nan(NumHit_act,1);
for iC = 1:NumHit_act  % loop through Hit-activated cholinergic cells
    inx = Hit_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    
    % Foreperiod, reward size, number of evoked spikes
    TE = loadcb(cellid,'TrialEvents');  % trial events
    goITI = TE.ITIDistribution;   % foreperiod
    goITI = goITI(~isnan(TE.Hit));
    goTITI = TE.TotalITI;   % realized foreperiod (with restarts)
    goTITI = goTITI(~isnan(TE.Hit));
    goRew = TE.LeftWaterVolume;   % reward size
    goRew = goRew(~isnan(TE.Hit));
    gospno = Hit_spno{inx};   % number of spikes fired in response to reinforcement
    
    % Plot association between foreperiod and number of evoked spikes and return p-value
    p = plotassoc(gospno,goITI);  % plot association between foreperiod and spike number and return p-value
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('Foreperiod')
    xlabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_ITI_VS_GOSPNO'];
        saveas(gcf,fnm)   % save bar graph
    end
    
    % Plot association between full foreperiod and number of evoked spikes and return p-value
    p = plotassoc(gospno,goTITI);  % plot association between full foreperiod and spike number and return p-value
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('Full foreperiod')
    xlabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_TITI_VS_GOSPNO'];
        saveas(gcf,fnm)   % save bar graph
    end
    
    % Plot association between foreperiod and number of evoked spikes and return p-value
    p = plotassoc(goRew,gospno);  % plot association between reward size and spike number and return p-value
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    xlabel('Reward size')
    ylabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_REW_VS_GOSPNO2'];
        saveas(gcf,fnm)   % save bar graph
    end
    
    goallp(inx) = p;  % p-values for regression
    close all
    waitbar(iC/NumHit_act)
end
close(wb)

% -------------------------------------------------------------------------
function p = plotassoc(X,Y)

% Association between X and Y
unq = sort(unique(X),'ascend');   % different values of X
NumSI = length(unq);   % number of different values
[xk yk] = deal([]);
[ngmn ngse] = deal(nan(1,NumSI));
for k = 1:NumSI  % loop through different values of X
    sic = unq(k);   % current value
    ngss = Y(X==sic);   % number of Y values corresponding to the given X
    nm = sum(X==sic);
    xk = [xk ngss];   % Y
    yk = [yk ones(1,nm)*k];   % grouping variable
    ngmn(k) = nanmean(ngss);   % mean
    ngse(k) = nanstd(ngss) / sqrt(sum(~isnan(ngss)));   % SE
end

% Plot
figure
subplot(1,2,1)
boxplot(xk,yk)
subplot(1,2,2)
bar(unq,ngmn,0.5,'FaceColor','none')
hold on
errorbar(unq,ngmn,ngse,'k+')

% Regression
y = Y;
x = X;
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
disp(p)