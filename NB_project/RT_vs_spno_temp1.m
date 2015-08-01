%%

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

% Association between RT/stimulus intensity and cholinergic response
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
    inx1 = ngspno==0;   % no response
    inx2 = ngspno>=1;   % response
%     if sum(inx1) * sum(inx2) > 0   % neither of the index sets are empty
%         boxstat(ngRT(inx1),ngRT(inx2),'no response','response');  % stat., false alarms
%         title([regexprep(cellid,'_',' ') ' FA'])
%     end
    
    % Number of eveoked spikes vs. stimulus intensity
    unq = sort(unique(ngspno),'ascend');   % different spike numbers of evoked bursts
    NumSpno = length(unq);   % number of different burst lengths
    [xk yk] = deal([]);
    [ngmn ngse] = deal(nan(1,NumSpno));
    for k = 1:NumSpno  % loop different spike numbers
        sic = unq(k);   % number of spikes in the evoked burst
        ngsis = ngSI(ngspno==sic);   % stim. intensity corresponding to the selected bursts
        nm = sum(ngspno==sic);
        xk = [xk ngsis];   % stim. intensities
        yk = [yk ones(1,nm)*k];   % grouping variable
        ngmn(k) = mean(ngsis);   % mean stim intensity
        ngse(k) = std(ngsis) / sqrt(nm);   % SE of stim. intensity
    end
    
    % Plot
    figure
    boxplot(xk,yk)
    figure
    bar(unq,ngmn,0.5,'FaceColor','none')
    hold on
    errorbar(unq,ngmn,ngse,'k+')
    
    % Regression
    y = ngspno;   % number of evoked spikes
    x = ngSI;   % stimulus intensity
    [b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
    R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
    F = stats(2);           % F-test for H0: all coeff.-s are zero
    p = stats(3);           % F-test significance
    disp(p)
    ngallp(inx) = p;  % p-values for regression
    keyboard
end

goallp = nan(NumHit_act,1);
for iC = 1:NumHit_act  % loop through Hit-activated cholinergic cells
    inx = Hit_act_inx(iC);
    goRT = TE.GoRT;   % go reaction time
    goRT = goRT(~isnan(goRT));
    goSI = TE.SoundIntensity;   % stim. intensity
    goSI = goSI(~isnan(TE.Hit));
    gospno = Hit_spno{inx};   % number of spikes fired in response to reinforcement
    inx1 = gospno==0;   % no response
    inx2 = gospno>=1;   % response
%     if sum(inx1) * sum(inx2) > 0   % neither of the index sets are empty
%         boxstat(goRT(inx1),goRT(inx2),'no response','response');  % stat., hits
%         title([regexprep(cellid,'_',' ') ' HIT'])
%     end
    
    % Regression
    y = gospno;
    x = goSI;
    [b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
    R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
    F = stats(2);           % F-test for H0: all coeff.-s are zero
    p = stats(3);           % F-test significance
%     disp(p)
    goallp(inx) = p;
    keyboard
end

%%

