%%

global FTHETA
ftheta = FTHETA;
global FNOTH
fnoth = FNOTH;
global DATINX1_THETA
datinx1_theta = DATINX1_THETA;
global DATINX2_THETA
datinx2_theta = DATINX2_THETA;
global DATINX1_NOTH
datinx1_noth = DATINX1_NOTH;
global DATINX2_NOTH
datinx2_noth = DATINX2_NOTH;

sr = 10000;
dsr = 1000;
const = sr / dsr;

% Plot
x1 = 0;
x2 = (maxi * winlen * const) / sr;
xa = linspace(x1,x2,length(aHx));
lna = floor((datinx2_theta-datinx1_theta)/winlen) * winlen / sr * 1000

%%

x = sum(data2(:,1:lna));
y = sum(data1(:,1:lna));
MC = 1:50:max(x);

next = 1;
fcn = [];
fsem = [];
fsd = [];
for mc = MC(2:end);
    xinx = find(x>MC(next)&x<MC(next+1));
    yc = y(xinx);
    fcn(next) = mean(yc)
    fsem(next) = std(yc) / sqrt(length(xinx));
    fsd(next) = std(yc);
    next = next + 1;
end

figure
plot((MC(2:end)+MC(1:end-1))/2,fcn)
hold on
errorbar((MC(2:end)+MC(1:end-1))/2,fcn,fsd,'k+')

%%

x = sum(data2(:,1:lna));
y = sum(data1(:,1:lna));
MC = 1:10:max(x);
winlen = 50;

next = 1;
fcn = [];
fsem = [];
fsd = [];
for mc = MC(5:end);
    xinx = find(x>MC(next)-winlen+1&x<MC(next));
    yc = y(xinx);
    fcn(next) = mean(yc)
    fsem(next) = std(yc) / sqrt(length(xinx));
    fsd(next) = std(yc);
    next = next + 1;
end

figure
plot((MC(5:end)+MC(1:end-4))/2,fcn)
hold on
errorbar((MC(5:end)+MC(1:end-4))/2,fcn,fsem,'k+')