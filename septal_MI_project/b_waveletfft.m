% Imported from Viktor.

% function b_waveletfft(wave)

clear all;
in2;
disc_futtat;
resamp=input('Downsampling frequency: ');
dj=input('Wavelet resolution (default is 0.05): ');
s0=2*(1/resamp);
step=fix(10000/resamp);
eeg=eeg(1:step:end);
dt=1/resamp;

dat=eeg;
variance=std(dat)^2;
dat=(dat-mean(dat))/sqrt(variance);
n=length(dat);
time=[0:length(dat)-1]*dt; 
pad=1;
j1=ceil((log(n*dt/s0)/log(2))/dj);
j=(0:j1);
s=s0.*2.^(j*dj);
omega0=6;
c=4*pi/(omega0+sqrt(2+omega0^2));
fperiod=c.*s;   
f=1./fperiod;  
lag1=0.72;  
mother='Morlet';
disp(' COMPUTING THE EEG WAVELET SPECTRUM ...')
[wave,period,scale,coi]=wavelet(dat,dt,pad,dj,s0,j1,mother);
wave_eeg=wave(:,1:end-1);
we=abs(wave_eeg).^2;

% figure
% sw1 = size(wave,1);
% subplot(3,1,1)
% w1 = wave(1:20,:);
% p1 = zeros(4097,20);
% f1 = zeros(4097,1);
% for i = 1:20
%     [p1(:,i) f1(:,1)] = pwelch(w1(i,:),[],[],8192,5000);
% end
% p1 = p1';
% f1 = f1';
% pet1 = p1(:,4:200);
% contourf_d(pet1);
% subplot(3,1,2)
% w2 = wave(21:40,:);
% p2 = zeros(4097,20);
% f2 = zeros(4097,1);
% for i = 1:20
%     [p2(:,i) f2(:,1)] = pwelch(w2(i,:),[],[],8192,5000);
% end
% p2 = p2';
% f2 = f2';
% pet2 = p2(:,4:200);
% contourf_d(pet2);
% subplot(3,1,3)
% w3 = wave(41:end,:);
% p3 = zeros(4097,20);
% f3 = zeros(4097,1);
% for i = 1:sw1-40
%     [p3(:,i) f3(:,1)] = pwelch(w3(i,:),[],[],8192,5000);
% end
% p3 = p3';
% f3 = f3';
% pet3 = p3(:,4:200);
% contourf_d(pet3);

figure
mx = max(we);
sw2 = size(we,2);
z = zeros(1,sw2);
for i = 1:sw2
    z(i) = find(we(:,i)==mx(i));
end
plot([1:sw2],z,'.')
set(gca,'ylim',[1 length(f)]);