function msptwocells2
%MSPTWOCELLS2   Analysis of MS cell pairs.
%   MSPTWOCELLS2 calculates auto-, crosscorrelation, and autospectrum for MS
%   cell pairs with theta-rhythmic firing.
%
%   See also MSPACORR.

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'MSHCsp\Viktor5\'];
resdir = [DATAPATH 'MSHCsp\Viktor5_2_43_45_46_22\'];
thetadir = [DATAPATH 'MSHCsp\Wavelet\theta_segments\'];
nonthetadir = [DATAPATH 'MSHCsp\Wavelet\nontheta_segments\'];
mm = pwd;

% Filelist
rid = '1237';
flist = {'201003291' '201003292' '201003293' '201003294'};
sf = length(flist);

% Main
sr = 20000;
o = 2;
fname = flist{o}     % 'filename'
ff = [thetadir 'THETA_SEGMENTS_' rid '_' fname '.mat'];     % load theta segments
load(ff)
ThetaSegments = uniteseg(ThetaSegments,sr);     % drop gaps < 0.5 s
ThetaSegments = short_killer(ThetaSegments);    % drop segments < 3 s
lent = ThetaSegments(2,:) - ThetaSegments(1,:);
minx = find(lent==max(lent));
th1 = ThetaSegments(1,minx);
th2 = ThetaSegments(2,minx);
    
ff = [nonthetadir 'NONTHETA_SEGMENTS_' rid '_' fname '.mat'];     % load non-theta segments
load(ff)
NonThetaSegments = short_killer(NonThetaSegments);    % drop segments < 3 s
lenn = NonThetaSegments(2,:) - NonThetaSegments(1,:);
minx = find(lenn==max(lenn));
nth1 = NonThetaSegments(1,minx);
nth2 = NonThetaSegments(2,minx);

shankno = 4;
ff = [inpdir fname '.res.' num2str(shankno)];       % load unit
res = textread(ff);     % spike timestamps
ff = [inpdir fname '.clu.' num2str(shankno)];
clu = textread(ff);     % indices of sorted cells
clu = clu(2:end);
ncells = max(clu);        % number of sorted cells (including multiunit)

% Autocorrelation
nc = 3;
vdisc = res(clu==nc)';
vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
if ~isempty(vdisc_theta)
    ac_theta = racorr(vdisc_theta/sr,500);    % theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_theta)),ac_theta)
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_THACORR.fig'];
    saveas(H,fns)
end
vdisc_nontheta = vdisc(vdisc>nth1&vdisc<nth2) - nth1;
if ~isempty(vdisc_nontheta)
    ac_nontheta = racorr(vdisc_nontheta/sr,500);    % non-theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_nontheta)),ac_nontheta)
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_NTHACORR.fig'];
    saveas(H,fns)
end
z1_theta = zeros(1,vdisc_theta(end)+10);
z1_theta(vdisc_theta) = 1;
z1_nontheta = zeros(1,vdisc_nontheta(end)+10);
z1_nontheta(vdisc_nontheta) = 1;
ac1_theta = ac_theta;
ac1_nontheta = ac_nontheta;
vdisc_theta1 = vdisc_theta;
vdisc_nontheta1 = vdisc_nontheta;

nc = 5;
vdisc = res(clu==nc)';
vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
if ~isempty(vdisc_theta)
    ac_theta = racorr(vdisc_theta/sr,500);    % theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_theta)),ac_theta,'FaceColor','red')
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_THACORR.fig'];
    saveas(H,fns)
end
vdisc_nontheta = vdisc(vdisc>nth1&vdisc<nth2) - nth1;
if ~isempty(vdisc_nontheta)
    ac_nontheta = racorr(vdisc_nontheta/sr,500);    % non-theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_nontheta)),ac_nontheta,'FaceColor','red')
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_NTHACORR.fig'];
    saveas(H,fns)
end
z2_theta = zeros(1,vdisc_theta(end)+10);
z2_theta(vdisc_theta) = 1;
z2_nontheta = zeros(1,vdisc_nontheta(end)+10);
z2_nontheta(vdisc_nontheta) = 1;
ac2_theta = ac_theta;
ac2_nontheta = ac_nontheta;
vdisc_theta2 = vdisc_theta;
vdisc_nontheta2 = vdisc_nontheta;

nc = 6;
vdisc = res(clu==nc)';
vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
if ~isempty(vdisc_theta)
    ac_theta = racorr(vdisc_theta/sr,500);    % theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_theta)),ac_theta,'FaceColor','green')
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_THACORR.fig'];
    saveas(H,fns)
end
vdisc_nontheta = vdisc(vdisc>nth1&vdisc<nth2) - nth1;
if ~isempty(vdisc_nontheta)
    ac_nontheta = racorr(vdisc_nontheta/sr,500);    % non-theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_nontheta)),ac_nontheta,'FaceColor','green')
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_NTHACORR.fig'];
    saveas(H,fns)
end
z3_theta = zeros(1,vdisc_theta(end)+10);
z3_theta(vdisc_theta) = 1;
z3_nontheta = zeros(1,vdisc_nontheta(end)+10);
z3_nontheta(vdisc_nontheta) = 1;
ac3_theta = ac_theta;
ac3_nontheta = ac_nontheta;
vdisc_theta3 = vdisc_theta;
vdisc_nontheta3 = vdisc_nontheta;

shankno = 2;
ff = [inpdir fname '.res.' num2str(shankno)];       % load unit
res = textread(ff);     % spike timestamps
ff = [inpdir fname '.clu.' num2str(shankno)];
clu = textread(ff);     % indices of sorted cells
clu = clu(2:end);
ncells = max(clu);        % number of sorted cells (including multiunit)
nc = 2;
vdisc = res(clu==nc)';
vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
if ~isempty(vdisc_theta)
    ac_theta = racorr(vdisc_theta/sr,500);    % theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_theta)),ac_theta,'FaceColor','magenta')
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_THACORR.fig'];
    saveas(H,fns)
end
vdisc_nontheta = vdisc(vdisc>nth1&vdisc<nth2) - nth1;
if ~isempty(vdisc_nontheta)
    ac_nontheta = racorr(vdisc_nontheta/sr,500);    % non-theta autocorrelation
    H = figure;
    bar(linspace(-1000,1000,length(ac_nontheta)),ac_nontheta,'FaceColor','magenta')
    fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_NTHACORR.fig'];
    saveas(H,fns)
end
z4_theta = zeros(1,vdisc_theta(end)+10);
z4_theta(vdisc_theta) = 1;
z4_nontheta = zeros(1,vdisc_nontheta(end)+10);
z4_nontheta(vdisc_nontheta) = 1;
ac4_theta = ac_theta;
ac4_nontheta = ac_nontheta;
vdisc_theta4 = vdisc_theta;
vdisc_nontheta4 = vdisc_nontheta;

% Raw data
H = figure;
plot(z2_nontheta,'r')
ylim([-2 2])
hold on
plot(z1_nontheta+1)
plot(z3_nontheta-1,'g')
plot(z4_nontheta-2,'m')
fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_RAW_NTH.fig'];
saveas(H,fns)

H = figure;
plot(z2_theta,'r')
ylim([-2 2])
hold on
plot(z1_theta+1)
plot(z3_theta-1,'g')
plot(z4_theta-2,'m')
fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_RAW_TH.fig'];
saveas(H,fns)

% Autospectrum
[y,w] = b_fft2(ac1_nontheta(1:500),200);    % non-theta
H = figure;
plot(w(2:25),y(2:25)/sum(y))
[y,w] = b_fft2(ac2_nontheta(1:500),200);
hold on
plot(w(2:25),y(2:25)/sum(y),'r')
[y,w] = b_fft2(ac3_nontheta(1:500),200);
plot(w(2:25),y(2:25)/sum(y),'g')
[y,w] = b_fft2(ac4_nontheta(1:500),200);
plot(w(2:25),y(2:25)/sum(y),'m')
fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_AUTO_NTH.fig'];
saveas(H,fns)

[y,w] = b_fft2(ac1_nontheta(1:500),200);
H = figure;
plot(w(5:25),y(5:25)/sum(y(5:25)))
[y,w] = b_fft2(ac2_nontheta(1:500),200);
hold on
plot(w(5:25),y(5:25)/sum(y(5:25)),'r')
[y,w] = b_fft2(ac3_nontheta(1:500),200);
plot(w(5:25),y(5:25)/sum(y(5:25)),'g')
[y,w] = b_fft2(ac4_nontheta(1:500),200);
plot(w(5:25),y(5:25)/sum(y(5:25)),'m')
fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_AUTO_NTH_NORM2.fig'];
saveas(H,fns)


[y,w] = b_fft2(ac1_theta(1:500),200);    % theta
H = figure;
plot(w(2:25),y(2:25)/sum(y))
[y,w] = b_fft2(ac2_theta(1:500),200);
hold on
plot(w(2:25),y(2:25)/sum(y),'r')
[y,w] = b_fft2(ac3_theta(1:500),200);
plot(w(2:25),y(2:25)/sum(y),'g')
[y,w] = b_fft2(ac4_theta(1:500),200);
plot(w(2:25),y(2:25)/sum(y),'m')
fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_AUTO_TH.fig'];
saveas(H,fns)

[y,w] = b_fft2(ac1_theta(1:500),200);
H = figure;
plot(w(5:25),y(5:25)/sum(y(5:25)))
[y,w] = b_fft2(ac2_theta(1:500),200);
hold on
plot(w(5:25),y(5:25)/sum(y(5:25)),'r')
[y,w] = b_fft2(ac3_theta(1:500),200);
plot(w(5:25),y(5:25)/sum(y(5:25)),'g')
[y,w] = b_fft2(ac4_theta(1:500),200);
plot(w(5:25),y(5:25)/sum(y(5:25)),'m')
fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc) '_AUTO_TH_NORM2.fig'];
saveas(H,fns)

% Crosscorrelation
for k1 = 1:4
    for k2 = (k1+1):4
        str = ['[H1 H2 trsc] = czxcorr(vdisc_nontheta' num2str(k1)...
            '''/20000,vdisc_nontheta' num2str(k2) '''/20000);'];
        eval(str)
        figure(H1)
        title([num2str(k1) '-' num2str(k2) ' nontheta'])
        fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc)...
            '_CROSS_NTH_' num2str(k1) num2str(k2) '.fig'];
        saveas(H1,fns)
        figure(H2)
        title([num2str(k1) '-' num2str(k2) ' nontheta'])
        fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc)...
            '_NCROSS_NTH_' num2str(k1) num2str(k2) '.fig'];
        saveas(H2,fns)
        
        str = ['[H1 H2 trsc] = czxcorr(vdisc_theta' num2str(k1)...
            '''/20000,vdisc_theta' num2str(k2) '''/20000);'];
        eval(str)
        figure(H1)
        title([num2str(k1) '-' num2str(k2) ' theta'])
        fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc)...
            '_CROSS_TH_' num2str(k1) num2str(k2) '.fig'];
        saveas(H1,fns)
        figure(H2)
        title([num2str(k1) '-' num2str(k2) ' theta'])
        fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc)...
            '_NCROSS_TH_' num2str(k1) num2str(k2) '.fig'];
        saveas(H2,fns)

        str = ['H = lxcorr(vdisc_nontheta' num2str(k1)...
            '''/20000,vdisc_nontheta' num2str(k2) '''/20000);'];
        eval(str)
        title([num2str(k1) '-' num2str(k2) ' nontheta'])
        fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc)...
            '_CROSS200_NTH_' num2str(k1) num2str(k2) '.fig'];
        saveas(H,fns)
        
        str = ['H = lxcorr(vdisc_theta' num2str(k1)...
            '''/20000,vdisc_theta' num2str(k2) '''/20000);'];
        eval(str)
        title([num2str(k1) '-' num2str(k2) ' theta'])
        fns = [resdir rid '_' fname '_' num2str(shankno) '_' num2str(nc)...
            '_CROSS200_TH_' num2str(k1) num2str(k2) '.fig'];
        saveas(H,fns)
    end
end

% H1 = lxcorr(vdisc_theta1'/20000,(vdisc_theta1'+2000)/20000);  % CCG test

1



% -------------------------------------------------------------------------
function segments2 = uniteseg(segments,sr)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    if (la-fi)/sr < 0.5
        [fnx fny] = find(segments2==fi);
        segments2(fnx,fny) = segments(2,k+1);
        segments2 = [segments2(1:2,1:fny) segments2(1:2,fny+2:end)];
    end
end

% ----------------------------------------------------------------------------------
function segments = short_killer(segments)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<30000);         % leaving segments shorter than 3 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];

% -------------------------------------------------------------------------
function sacr = racorr(ncc,bno)
%RACORR   Autocorrelation.
%   RACORR2(VD,BNO) calculates autocorrelogram for discriminated unit VD, 
%   using a +-1000 ms time window and BNO bins.

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;

% Autocorrelogram
zunit1 = zeros(1,length(round(nc))+5);
zunit1(round(nc)) = 1;
acr = xcorr(zunit1,1*sr);
acr(length(acr)/2+0.5) = [];
acr = reshape(acr,length(acr)/bno,bno);     % window: -200 ms - 200 ms
sacr = sum(acr);

% -------------------------------------------------------------------------
function H1 = lxcorr(ncc1,ncc2)
%LXCORR   Crosscorrelation.
%   LXCORR(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-200 ms time window.
%
%   H1 = LXCORR(VD1,VD2) returns the handle of the resulting plot.
%
%   See also CZXCORR, XCORR and CZACORR.

% Input argument check
error(nargchk(2,2,nargin))

% Calculate spike times in milliseconds
sr = 1000;
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
nc1(nc1<0.5) = [];
nc2(nc2<0.5) = [];

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
ccr = xcorr(zunit2,zunit1,0.2*sr);     % 1->2; window: -50 ms - 50 ms
ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
H1 = figure;
bar(linspace(-200,200,length(ccr)),ccr)
set(gca,'XLim',[-200 200])