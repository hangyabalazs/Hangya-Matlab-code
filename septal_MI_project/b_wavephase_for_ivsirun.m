function [ds4,t3] = b_wavephase_for_ivsirun
%WAVEPHASE_FOR_IVSIRUN Verson of wavephase used by IVSI.
%   WAVEPHASE_FOR_IVSIRUN has two output arguments: one contains the localization of the first
%   and last points of the theta intervals of the eeg, the other contains the localization 
%   of the eeg minimums.
%
%   WAVEPHASE_FOR_IIVSIRUN uses the following theta definition: more than 30 % of the eeg
%   wavelet power fall in the 3 - 6 Hz band.
%
%   Imported from Viktor.
%
%   See also WAVEPHASE, IVSI and IVSIRUN.

%Input arguments check
error(nargchk(0,0,nargin));

%Input variables
global IN
data = IN{1};
eeg = IN{2};
fname = IN{3};
pathname = IN{4};
datinx1 = IN{5};
datinx2 = IN{6};
time = IN{7};
unit = IN{8};
dt = IN{9};
meret = IN{10};
mintafr = IN{11};
xlimit = IN{12};

%Wavelet spectrum calculation
dat5=unit;
leneeg=length(eeg); 
dat=eeg(1:20:leneeg); 
variance=std(dat)^2;
dat=(dat-mean(dat))/sqrt(variance) ;
n=length(dat);
dt=0.02;							% for sample frequency of 500 Hz
dtunit=0.0001;                      % for sample frequency of 10000 Hz
time=[0:length(dat)-1]*dt; 
xlim=[min(time),max(time)];  
timeunit=[0:length(dat5)-1]*dtunit;
xlimunit=[min(timeunit),max(timeunit)];
pad=1;      
dj=0.1;    
s0=1*dt;    
j1=10/dj;    
lag1=0.72;  
mother='Morlet';
disp(' COMPUTING THE WAVELET SPECTRUM ...')
[wave,period,scale,coi]=wavelet(dat,dt,pad,dj,s0,j1,mother);
power=(abs(wave)).^2 ;        
f1=1; f2=30;
pf1=1/f1; pf2=1/f2;
i2=max(find(period<=pf1)); i1=max(find(period<=pf2));
ff1=2; ff2=8; ff3=10; ff4=20;
pff1=1/ff1; pff2=1/ff2; pff3=1/ff3; pff4=1/ff4;
if2=max(find(period<=pff1)); if1=max(find(period<=pff2));
if4=max(find(period<=pff3)); if3=max(find(period<=pff4));
pow=power(i1:i2,:); pow=pow';			 
globws=variance*(sum(pow)/n);   		 
pow1=power(if1:if2,:); pow1=pow1';	 	 
glwsstd=variance*std(pow);
totws=variance*mean(pow1);
totwsstd=variance*std(pow1);
pow2=power(if3:if4,:); pow2=pow2';	 	
totws2=variance*mean(pow2);
totwstd2=variance*std(pow2);
Cdelta=0.776;   							
avg=find((scale>=1/ff2)&(scale<1/ff1));
scaleavg=(scale')*(ones(1,n));  
scaleavg=power./scaleavg;   
scaleavg=variance*dj*dt/Cdelta*sum(scaleavg(avg,:));  
scav=mean(scaleavg);
scavstd=std(scaleavg);
avg2=find((scale>=1/ff4)&(scale<1/ff3));
scalavg2=(scale')*(ones(1,n));  
scalavg2=power./scalavg2;   
scalavg2=variance*dj*dt/Cdelta*sum(scalavg2(avg2,:));  
scav2=mean(scaleavg);
scavstd2=std(scaleavg);

%Selection of segments w. >50% ~3-6 Hz frequency content
thp=zeros(1,n);
allp=zeros(1,n);
selector=zeros(1,n);
for i=1:n
    thp(1,i)=sum(power(60:80,i));
    allp(1,i)=sum(power(:,i));
end;
thprop=thp./allp;
for i=1:n
    if thprop(1,i)>0.5
        selector(1,i)=1;
    end;
end;
ds=diff(selector);
ds1=zeros(1,n-1);
ds1=abs(ds);
if selector(1)==1
    ds1(1)=1;
end;
if selector(end)==1
    ds1(end)=1;
end;
ds2=find(ds1(1,:)==1);
if rem(length(ds2),2)==1,
    error('Length ds2 is not evan!');
end;

e1=diff(ds2);   %leaving the intervals shorter then 1,5 sec...
e2=e1(1:2:end);
e3=find(e2>15000/20);
e4=zeros(1,2*length(e3));
for i=1:length(e3),
    e4(2*i)=2*e3(i);
    e4(2*i-1)=e4(2*i)-1;
end;
e5=ds2(e4);
ds3=e5;

ds4=ds3.*20;    %localization of first and last points of theta intervals
ds5=ds4(1,1:2:end); %first points of theta intervals
ds6=ds4(1,2:2:end); %last points of theta intervals
if length(ds5)>length(ds6)
    ds5=ds5(1:length(ds5)-1);
end;
if length(ds6)>length(ds5)
    ds6=ds6(1:length(ds6)-1);
end;
ds7=ds6-ds5;

eegs=zeros(1,length(eeg));  %HB
for i=1:length(eegs),   %HB
    f=find(ds4<i);  %HB
    if rem(length(f),2)==1, %HB
        eegs(i)=eeg(i); %HB
    end;    %HB
end;    %HB

%ds8=sum(ds7);   %HB
%eegs=zeros(1,ds8);  %HB
%next=1; %HB
%for i=1:length(ds7),    %HB
%    for j=1:ds7(i), %HB
%        eegs(next)=eeg(ds5(i)+j);   %HB
%        next=next+1;    %HB
%    end;    %HB
%end;    %HB

if length(eegs)>6150 fo=2048;
else fo=fix(length(eegs)/3);
end;
b=fir1(fo,0.0014);
theta=filtfilt(b,1,eegs);
figure
plot(theta)
%disp([datinx1+bg datinx1+fin])

%Minimums in EEG signal
pret=zeros(1,length(theta));
avtheta=mean(theta);
for i=1:length(theta)
    if theta(1,i)<avtheta
        pret(1,i)=theta(1,i);
    end;
end;
t=[1:length(theta)];
t=[t;pret(t(1:length(t)))];
t(1,length(t)+1)=length(pret)+2;
t1=diff(t(2,:));
t2=find(t1(1,1:length(t1)-1)<0 & t1(1,2:length(t1))>0);
t3=t(1,t2+1);   %here are the minimums
time=[1:length(theta)];
figure
plot(theta,'m')
set(gca,'XLim',[0 length(theta)])
hold on;
mi=min(theta); ma=max(theta);
for j=1:length(t3)
    rajz=[time(t3(j)) time(t3(j));mi ma];
    line(rajz(1,:),rajz(2,:));
end;
hold off;
title('Filtered EEG, minimums marked');

%   List of files calling WAVEPHASE_FOR_IVSIRUN:
%   b_ivsi_main_for_ivsirun