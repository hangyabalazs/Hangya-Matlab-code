%   Imported from Viktor.

% display the cross-wavelet abs. & phase

disp('non-theta interval')
in2;
nt_start=datinx1;
nt_stop=datinx2;
w_nt_start=fix(nt_start/25);
w_nt_stop=fix(nt_stop/25);
b1=fir1(2048,0.0014);
eeg_nt=eeg;
unit_nt=unit;
f_eeg_nt=filtfilt(b1,1,eeg_nt);

disp('theta interval')
in2;
t_start=datinx1;
t_stop=datinx2;
w_t_start=fix(t_start/25);
w_t_stop=fix(t_stop/25);
b2=fir1(2048,0.0014);
eeg_t=eeg;
unit_t=unit;
f_eeg_t=filtfilt(b2,1,eeg_t);
figure; plot(f_eeg_nt)
set(gca,'XLim',[1 length(f_eeg_nt)]);
figure; plot(unit_nt)
set(gca,'XLim',[1 length(unit_nt)]);
figure; plot(f_eeg_t)
set(gca,'XLim',[1 length(f_eeg_t)]);
figure; plot(unit_t)
set(gca,'XLim',[1 length(unit_t)]);

disp('cross wavelets')

k=input('Which file? ');
cmnd=['load c:\_analysis\nora\cross_wavelet\kdat' num2str(k) 'cross'];
eval(cmnd)
clear cmnd 
w_segm=input('Which segment? ');
timecw=[(1/400):(1/400):(12000/400)];
n1=length(wavcross1);
dt = 0.0025 ;
dj = 0.05;    
s0 = 2*dt;    
j1 = fix(7/dj);
fourier_factor = (4*pi)/(6.0 + sqrt(2 + 6.0^2));
scale=zeros(1,j1);
frw=zeros(1,j1);
for i=1:j1+1
    scale(i)=fourier_factor *s0*2.^((i-1)*dj);
    frw(i)=1/scale(i);
end;
l_segm=w_nt_stop-(w_nt_start+1);
l_segm1=w_t_stop-(w_t_start+1);
levels=2.^[-4:0.2:7];
w_t_start=w_t_start-(600000/25); w_t_stop=w_t_stop-(600000/25);
wabs3=abs(wavcross3);

switch w_segm
case 1
    wabs1=abs(wavcross1);
    wabs1=flipud(wabs1(97:141,:));
    maxvalue_of_all = max(max([wabs1(:,w_nt_start:w_nt_stop) wabs3(97:141,w_t_start+1:w_t_stop)]));
    maxvalue_of_act = max(max(wabs1(:,w_nt_start:w_nt_stop)));
    rate = maxvalue_of_act / maxvalue_of_all;
    maxcolor = fix(rate * 256);
    C = b_defaultmap(8);
    C = C(1:maxcolor,:);
    figure;
    colormap(C);
    contourf_d(wabs1(:,w_nt_start:w_nt_stop),levels);
case 2
    wabs2=abs(wavcross2);
    wabs2=flipud(wabs2(97:141,:));
    w_nt_start=w_nt_start-(300000/25); w_nt_stop=w_nt_stop-(300000/25);
    maxvalue_of_all = max(max([wabs2(:,w_nt_start:w_nt_stop) wabs3(97:141,w_t_start+1:w_t_stop)]));
    maxvalue_of_act = max(max(wabs2(:,w_nt_start:w_nt_stop)));
    rate = maxvalue_of_act / maxvalue_of_all;
    maxcolor = fix(rate * 256);
    C = b_defaultmap(8);
    C = C(1:maxcolor,:);
    figure;
    colormap(C);
    contourf_d(wabs2(:,w_nt_start:w_nt_stop),levels);
case 3
    wabs3=abs(wavcross3);
    wabs3=flipud(wabs3(97:141,:));
    maxvalue_of_all = max(max([wabs3(:,w_t_start:w_t_stop) wabs3(97:141,w_t_start+1:w_t_stop)]));
    maxvalue_of_act = max(max(wabs3(:,w_t_start:w_t_stop)));
    rate = maxvalue_of_act / maxvalue_of_all;
    maxcolor = fix(rate * 256);
    C = b_defaultmap(8);
    C = C(1:maxcolor,:);
    figure;
    colormap(C);
    contourf_d(wabs3(:,w_t_start:w_t_stop),levels);
case 4
    wabs4=abs(wavcross4);
    wabs4=flipud(wabs4(97:141,:));
    w_nt_start=w_nt_start-(900000/25); w_nt_stop=w_nt_stop-(900000/25);
    maxvalue_of_all = max(max([wabs4(:,w_nt_start:w_nt_stop) wabs3(97:141,w_t_start+1:w_t_stop)]));
    maxvalue_of_act = max(max(wabs4(:,w_nt_start:w_nt_stop)));
    rate = maxvalue_of_act / maxvalue_of_all;
    maxcolor = fix(rate * 256);
    C = b_defaultmap(8);
    C = C(1:maxcolor,:);
    figure;
    colormap(C);
    contourf_d(wabs4(:,w_nt_start:w_nt_stop),levels);
case 5
    wabs5=abs(wavcross5);
    wabs5=flipud(wabs5(97:141,:));
    w_nt_start=w_nt_start-(1200000/25); w_nt_stop=w_nt_stop-(1200000/25);
    maxvalue_of_all = max(max([wabs5(:,w_nt_start:w_nt_stop) wabs3(97:141,w_t_start+1:w_t_stop)]));
    maxvalue_of_act = max(max(wabs5(:,w_nt_start:w_nt_stop)));
    rate = maxvalue_of_act / maxvalue_of_all;
    maxcolor = fix(rate * 256);
    C = b_defaultmap(8);
    C = C(1:maxcolor,:);
    figure;
    colormap(C);
    contourf_d(wabs5(:,w_nt_start:w_nt_stop),levels);
    %     case 6
    %         wabs6=abs(wavcross6);
    %         wabs6(97,w_nt_start)=maxcoeff;
    %         w_nt_start=w_nt_start-(1500000/25); w_nt_stop=w_nt_stop-(1500000/25);
    %         wabss6=[wabs6(97:141,w_nt_start:w_nt_stop) wabs3(97:141,w_t_start+1:w_t_stop)];
    %         wabss6=flipud(wabss6);
    %         figure; contourf_d(wabss6,levels);
end;

    wabs3=flipud(wabs3(97:141,:));
    maxvalue_of_all = max(max([wabs3(:,w_t_start:w_t_stop) wabs3(:,w_t_start+1:w_t_stop)]));
    maxvalue_of_act = max(max(wabs3(:,w_t_start:w_t_stop)));
    rate = maxvalue_of_act / maxvalue_of_all;
    maxcolor = fix(rate * 256);
    C = b_defaultmap(8);
    C = C(1:maxcolor,:);
    figure;
    colormap(C);
    contourf_d(wabs3(:,w_t_start:w_t_stop),levels);