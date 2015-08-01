% instfrek = [];lenu=length(eeg);
% isi = diff(vdisc);
% for i = 1:length(vdisc)-1
%     instfrek(vdisc(i):vdisc(i+1)) = 1 / isi(i);
% end
% instfrek(1:vdisc(1)-1) = 1 / vdisc(1);
% instfrek(vdisc(end):lenu) = 1 / (lenu - vdisc(end));
% instfrek2 = instfrek * 10000;
% figure;plot(instfrek2)

[powunit,phaseunit,f] = unitwavelet(vdisc,lenu);
figure
imagesc(powunit)
b_rescaleaxis('Y',f)
[poweeg,phaseeeg,f] = eegwavelet(eeg(1:10:end));
figure
imagesc(poweeg)
b_rescaleaxis('Y',f)

flt = fir1(4096,10/5000,'low');
feeg2 = filtfilt(flt,1,eeg);
flt = fir1(4096,50/5000,'low');
feeg = filtfilt(flt,1,eeg);

vd10 = vdisc(find(vdisc<300000&vdisc>200000));
isi10 = diff(vd10);
cnt10 = vd10(1:end-1) + round(isi10/2);
fcn = feeg(cnt10);
icn = instfrek(cnt10);
figure
plot(fcn,icn*10000,'.k')
icn = icn*10000;
[R,p] = corr(fcn',icn')
title('EEG vs. instfreq.')
xlabel('EEG')
ylabel('instfreq.')


T = [-25000:10:25000];
cc = zeros(1,length(T));
next = 1;
for t = T
    cnt2 = cnt10(find(cnt10+t>0));
    fcn2 = feeg2(cnt2+t);
    icn2 = instfrek(cnt2) * 10000;
    cc(next) = corr(fcn2',icn2');
    next = next + 1;
end
figure
plot(T,cc)
xlabel('time shift')
ylabel('instrfek-eeg corr')


isiloc = isi_vshort;    % or isi_sort, isi10
rk = tiedrank(isiloc);
nth = rk(1:end-1);
n1th = rk(2:end);
figure;
plot(nth,n1th,'.k')
xlabel('N. ISI')
ylabel('N+1. ISI')
title('Ranked recurrence plot')

nth = isiloc(1:end-1);
n1th = isiloc(2:end);
figure;
plot(nth,n1th,'.k')
xlabel('N. ISI')
ylabel('N+1. ISI')
title('Recurrence plot')


fn = find(isi>10000);
vd_short = vdisc;
vd_short(fn) = [];
isi_short = isi;
isi_short(fn) = [];
figure;
hist(isi_short,25)
fn2 = find(isi>2000);
isi_vshort = isi;
isi_vshort(fn2) = [];
figure;
hist(isi_vshort,25)

figure;
subplot(211)
plot(eeg(580000:650000))
subplot(212)
plot(unit(580000:650000))

figure;
subplot(211)
plot(eeg(1480000:1550000))
subplot(212)
plot(unit(1480000:1550000))

fis = find(isi>10000);
lenfis = length(fis);
for k = 1:lenfis
    if k == 1
        fireseg(1,k) = vdisc(1);
    else
        fireseg(1,k) = vdisc(fis(k-1)+1);
    end
    if k == length(fis)
        fireseg(2,k) = vdisc(end);
    else
        fireseg(2,k) = vdisc(fis(k));
    end
end

zs = zeros(1,lenfis);
os = ones(1,lenfis);
patcha = [fireseg(1,:); fireseg(1,:); fireseg(2,:); fireseg(2,:); fireseg(1,:)];
patchb = [zs; os; os; zs; zs];
figure
subplot(2,1,1)
ptch = patch(patcha,patchb,'red');
set(ptch,{'edgecolor'},get(ptch,{'facecolor'}));
subplot(2,1,2)
lent = size(ThetaSegments,2);
zs2 = zeros(1,lent);
os2 = ones(1,lent);
patcha2 = [ThetaSegments(1,:); ThetaSegments(1,:); ThetaSegments(2,:); ThetaSegments(2,:); ThetaSegments(1,:)];
patchb2 = [zs2; os2; os2; zs2; zs2];
ptch2 = patch(patcha2,patchb2,'blue');
set(ptch2,{'edgecolor'},get(ptch2,{'facecolor'}));

lent = size(ThetaSegments,2);
thcum = [];
for k = 1:lent
    thcum = [thcum ThetaSegments(1,k):ThetaSegments(2,k)];
end
frscum = [];
for k = 1:lenfis
    frscum = [frscum fireseg(1,k):fireseg(2,k)];
end
olp = length(intersect(thcum,frscum))

thonly = setdiff(thcum,frscum);
frsonly = setdiff(frscum,thcum);
thfrs = intersect(thcum,frscum);

ip = thonly;
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
prepa1 = zeros(2,lenfdr+1);
prepa1(1,1) = ip(1);
for t = 1:lenfdr
    prepa1(2,t) = ip(fdr(t));
    prepa1(1,t+1) = ip(fdr(t)+1);
end
prepa1(2,end) = ip(end);
zs = zeros(1,lenfdr+1);
os = ones(1,lenfdr+1);
patcha = [prepa1(1,:); prepa1(1,:); prepa1(2,:); prepa1(2,:); prepa1(1,:)];
patchb = [zs; os; os; zs; zs];
figure
ptch = patch(patcha,patchb,'blue');
set(ptch,{'edgecolor'},get(ptch,{'facecolor'}));

ip = frsonly;
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
prepa2 = zeros(2,lenfdr+1);
prepa2(1,1) = ip(1);
for t = 1:lenfdr
    prepa2(2,t) = ip(fdr(t));
    prepa2(1,t+1) = ip(fdr(t)+1);
end
prepa2(2,end) = ip(end);
zs = zeros(1,lenfdr+1);
os = ones(1,lenfdr+1);
patcha = [prepa2(1,:); prepa2(1,:); prepa2(2,:); prepa2(2,:); prepa2(1,:)];
patchb = [zs; os; os; zs; zs];
ptch = patch(patcha,patchb,'yellow');
set(ptch,{'edgecolor'},get(ptch,{'facecolor'}));

ip = thfrs;
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
prepa3 = zeros(2,lenfdr+1);
prepa3(1,1) = ip(1);
for t = 1:lenfdr
    prepa3(2,t) = ip(fdr(t));
    prepa3(1,t+1) = ip(fdr(t)+1);
end
prepa3(2,end) = ip(end);
zs = zeros(1,lenfdr+1);
os = ones(1,lenfdr+1);
patcha = [prepa3(1,:); prepa3(1,:); prepa3(2,:); prepa3(2,:); prepa3(1,:)];
patchb = [zs; os; os; zs; zs];
ptch = patch(patcha,patchb,'green');
set(ptch,{'edgecolor'},get(ptch,{'facecolor'}));

olp = length(thfrs);
tho = length(thonly);
frso = length(frsonly);
olp_pc = olp / lenu
tho_pc = tho / lenu
frso_pc = frso / lenu
non_pc = 1 - (olp_pc + tho_pc + frso_pc)

z = 0;
for k = 1:length(vdisc)
    fn1 = find(ThetaSegments<vdisc(k),1,'last');
    fn2 = find(ThetaSegments>vdisc(k),1,'first');
    if ~isempty(fn1)
        if isequal(rem(fn1,2),1)
            z = z + 1;
        end
    else
        if isequal(rem(fn2,2),0)
            z = z + 1;
        end
    end
end

flt = fir1(4096,[95/5000 145/5000]);
feeg_sharp = filtfilt(flt,1,eeg);
figure;plot(feeg_sharp)
