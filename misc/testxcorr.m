nc1 = randpoisson(5000,1000000);   % 5000 spikes in 1000 s
nc1(nc1>1000000) = [];
nc1(nc1<0.5) = [];

brst = [0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 zeros(1,35)];

zunit1 = zeros(1,1000000);
zunit1(round(nc1)) = 1;
zunit2 = repmat(brst,1,1000000/50);

wn = 1000;
ccr = xcorr(zunit2,zunit1,1000);

H1 = figure;
bar(linspace(-wn,wn,length(ccr)),ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])