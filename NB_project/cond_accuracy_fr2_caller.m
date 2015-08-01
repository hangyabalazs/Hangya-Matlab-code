%% call 'cond_accuracy_fr2'

figure
T = -60:1:-2;
magn = nan(1,length(T));
cnt = 0;
for lim1 = T
    cnt = cnt + 1;
    lim2 = lim1 + 2;
    condPerf = cond_accuracy_fr2(lim1,lim2);
    magn(cnt) = condPerf(3);
    hold on
end

figure
plot(T,magn)