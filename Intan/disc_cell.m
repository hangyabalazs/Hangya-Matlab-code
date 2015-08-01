%% plot

figure;plot(t(1:100:1000000),data(1:100:1000000,14))
figure;plot(t(1:10000),data(1:10000,13:16)+repmat(0:40:120,10000,1))

%% discriminate

D = data(1:10000,14);
vdisc = b_udisc(-D);