%%

for k = 1:length(wmax)
    df(k) = circdiff(orig_maxphase(k),wmax(k),'deg');
end

figure
plot(df,'.')

%%

mean(abs(df))
median(abs(df))
circular_mean(abs(df),'deg')