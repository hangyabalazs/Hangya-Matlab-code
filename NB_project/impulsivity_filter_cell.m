%% function impulsivity_filter

[Bl, Ia, Ib] = unique(TE.BlockNum);
bGoPerf = nan(1,length(Bl));
for ibl = 1:length(Bl),
    bGoPerf(ibl) = nansum(TE.Hit(Ib==ibl))/(nansum(TE.Hit(Ib==ibl))+nansum(TE.Miss(Ib==ibl)));
    bNoGoPerf(ibl) = nansum(TE.FalseAlarm(Ib==ibl))/(nansum(TE.FalseAlarm(Ib==ibl))+nansum(TE.CorrectRejection(Ib==ibl)));
end

figure
plot(bNoGoPerf)


%%
NUMtrials = length(TE.TrialStart);
win = 10;
for ibl = 1:NUMtrials-win
    bGoPerf(ibl) = nansum(TE.Hit(ibl:ibl+win))/(nansum(TE.Hit(ibl:ibl+win))+nansum(TE.Miss(ibl:ibl+win)));
    bNoGoPerf(ibl) = nansum(TE.FalseAlarm(ibl:ibl+win))/(nansum(TE.FalseAlarm(ibl:ibl+win))+nansum(TE.CorrectRejection(ibl:ibl+win)));
    bLengthITI(ibl) = mean(TE.TotalITI(ibl:ibl+win));
end

figure
plot(bNoGoPerf)

figure
plot(bLengthITI)

figure;plot(TE.FalseAlarm.*TE.StimulusDuration/10,'r.','MarkerSize',16)
hold on
plot(TE.CorrectRejection.*TE.StimulusDuration/10,'g.','MarkerSize',16)


%%

allicks = [];
ILI = cell(1,NUMtrials);
mILI = nan(1,NUMtrials);
for k = 1:NUMtrials
    allicks = [allicks; TE.LickIn{k}+TE.TrialStart(k)];
    ILI{k} = diff(TE.LickIn{k});
    mILI(k) = mean(ILI{k});
end

figure
plot(mILI)