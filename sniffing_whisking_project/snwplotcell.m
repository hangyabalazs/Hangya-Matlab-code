%%

load('C:\Balazs\_data\HSW\R7\2008-09-03_16-51-55\SIGNALS.mat')

%%

figure;plot(-H1)
hold on
plot(-SN*10,'r')
plot(-W1,'g')

%%

load('C:\Balazs\_data\HSW\R3\2008-07-09_14-46-47_R3\ALIGNED_SIGNALS.mat')

%%

PELLET_H1 = PELLET_H1';
PELLET_SN = PELLET_SN';
PELLET_W1 = PELLET_W1';
figure;plot(-PELLET_H1(:))
hold on
plot(PELLET_SN(:),'r')
plot(-PELLET_W1(:),'g')