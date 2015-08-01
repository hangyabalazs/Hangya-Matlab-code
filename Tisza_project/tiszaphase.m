function tiszaphase

% Directories
inpdir = 'd:\_analysis\matlab_data\Tisza\';
resdir = 'd:\_analysis\matlab_data\Tisza\';

% Import
load([inpdir 'tivadq_interp.mat'])
x_tivad = xi;
y_tivad = yi;
ynorm_tivad = (y_tivad - mean(y_tivad)) / std(y_tivad);
load([inpdir 'csengerq_interp.mat'])
x_csenger = xi;
y_csenger = yi;
ynorm_csenger = (y_csenger - mean(y_csenger)) / std(y_csenger);

% Crosswavelet
sfr = 96;       % sampling frequency: 96/day
[wave_tiv,f] = wavelet(y_tivad,sfr);
[wave_csen,f] = wavelet(y_csenger,sfr);
wave_cross = wave_tiv .* conj(wave_csen);
figure
imagesc(abs(wave_tiv).^2)
b_rescaleaxis('Y',f)
saveas(gca,[resdir 'power_tiv.jpg'])
imagesc(abs(wave_csen).^2)
b_rescaleaxis('Y',f)
saveas(gca,[resdir 'power_csen.jpg'])
imagesc(abs(wave_cross).^2)
b_rescaleaxis('Y',f)
saveas(gca,[resdir 'power_cross.jpg'])
imagesc(angle(wave_tiv))
% b_rescaleaxis('Y',f)
saveas(gca,[resdir 'angle_tiv.jpg'])
imagesc(angle(wave_csen))
% b_rescaleaxis('Y',f)
saveas(gca,[resdir 'angle_csen.jpg'])
imagesc(angle(wave_cross))
% b_rescaleaxis('Y',f)
saveas(gca,[resdir 'angle_cross.jpg'])
clear wave_tiv wave_csen

% Discrimination
y_tiv = (y_tivad-mean(y_tivad)) / std(y_tivad);
y_csen = (y_csenger-mean(y_csenger)) / std(y_csenger);
flt = fir1(1024,(12/365)/(sfr/2),'high');
fy_tiv = filtfilt(flt,1,y_tiv);
fy_csen = filtfilt(flt,1,y_csen);
thres1 = 0.5;
vdtiv = discrim(fy_tiv,thres1);
thres2 = 0.5;
vdcsen = discrim(fy_csen,thres2);

% Crosswavelet angles
limit1 = 69;
limit2 = 85;
limit3 = 95;
limit4 = 106;
[cang_tc1 cmvl_tc1 angvec_tc1] = cangle(vdtiv,wave_cross,limit1,limit2,f);
[cang_ct1 cmvl_ct1 angvec_ct1] = cangle(vdcsen,wave_cross,limit1,limit2,f);
[cang_tc2 cmvl_tc2 angvec_tc2] = cangle(vdtiv,wave_cross,limit2,limit3,f);
[cang_ct2 cmvl_ct2 angvec_ct2] = cangle(vdcsen,wave_cross,limit2,limit3,f);
[cang_tc3 cmvl_tc3 angvec_tc3] = cangle(vdtiv,wave_cross,limit3,limit4,f);
[cang_ct3 cmvl_ct3 angvec_ct3] = cangle(vdcsen,wave_cross,limit3,limit4,f);
[L_tc2, U_tc2] = wconf(angvec_tc2,cang_tc2,cmvl_tc2)
[L_tc3, U_tc3] = wconf(angvec_tc3,cang_tc3,cmvl_tc3)



% -------------------------------------------------------------------------
function [wave,f] = wavelet(eeg,sfr)

% Filtering
% flt = fir1(512,[0.0027/(sfr/2) 24/(sfr/2)]);        % 1 hour - 1 year period time
% feeg = filtfilt(flt,1,eeg);
feeg = eeg;

% Standardization
variance = std(feeg) ^ 2;
feeg = (feeg - mean(feeg)) / sqrt(variance) ;

% Prepare for wavelet transformation
n = length(feeg);
dt = 1 / sfr;
pad = 1;
dj = 0.1;
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mis = -1;
% mif = 0.5;          %minimal intresting frequency
% mis = find(f>mif);
% mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(feeg,dt,pad,dj,s0,j1,mother,param,mis);



% -------------------------------------------------------------------------
function vdisc = discrim(unit,kuszob)

disc = find(unit>=kuszob);
discl = length(disc);
disc = [disc; unit(disc(1:discl))];
disc(1,discl+1) = length(unit) + 2;
dif = diff(disc(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j = 1:length(difn1) - 1
    [maxe,maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end



% -------------------------------------------------------------------------
function [cang cmvl wphh] = cangle(vdisc,wave_cross,pwind1,pwind2,f)

% Phase angles - crosswavelet
wph = angle(wave_cross(pwind1:pwind2,:));
wabs = abs(wave_cross(pwind1:pwind2,:)).^2;
mwabs = max(wabs);
lenw = length(wabs);

wph0 = zeros(1,lenw);
for k = 1:lenw
    mloc = find(wabs(:,k)==mwabs(k));
    wph0(k) = wph(mloc,k);
end
wphh = wph0(vdisc);
wahh = mwabs(vdisc);

% dwphh = diff(wphh);    % get "different" phase values
% fdw = [0 abs(dwphh)>0.05];
% wnew = [];
% for i = 1:length(fdw)
%     switch fdw(i)
%         case 0
%             if exist('act')
%                 act(end+1) = wphh(i);
%             else
%                 act = wphh(i);
%             end
%         case 1
%             if i>1 && fdw(i-1)==0
%                 wnew(end+1) = cmean(act);
%                 act = [];
%             elseif i>1 && fdw(i-1)==1
%                 wnew(end+1) = wphh(i-1);
%             end
%     end
% end
% wphh = wnew;

threshold = 10;     % drop subthreshold values
wphh = wphh(find(wahh>threshold));

n = length(wphh);
ftm = sum(exp(1).^(j*wphh)) / n;    % first trigonometric moment
cang = angle(ftm);   % mean angle
cmvl = abs(ftm);     % mean resultant length

cang / (2 * pi) / f(pwind1)
cang / (2 * pi) / f(pwind2)
cang / (2 * pi) / f(round((pwind1+pwind2)/2))



% -------------------------------------------------------------------------
function M = cmean(A)

ftm = sum(exp(1).^(j*A)) / length(A);    % first trigonometric moment
M = angle(ftm);   % mean angle