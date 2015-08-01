%% Import

n1=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_0.txt');
n2=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_1.txt');
n3=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_2.txt');
n4=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_3.txt');
n5=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_4.txt');
n6=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_5.txt');
n7=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_6.txt');
n8=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_7.txt');
n9=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\nr_01_8.txt');
ntest=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\Ntest_05b.txt');

%% Read position data

x_pos=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\x1_pos.txt')';
x_pos_ts=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\x1_pos_ts.txt')';
y_pos=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\y1_pos.txt')';
y_pos_ts=textread('F:\balazs\_analysis\Czurko\discriminated\acin06s022_sc2_Kall2\y1_pos_ts.txt')';

%% Sampling rate: 1000 Hz

sr = 1000;
n1_ = n1 * sr;    % in ms
n2_ = n2 * sr;
n3_ = n3 * sr;
n4_ = n4 * sr;
n5_ = n5 * sr;
n6_ = n6 * sr;
n7_ = n7 * sr;
n8_ = n8 * sr;
n9_ = n9 * sr;
ntest_ = ntest * sr;


%% Auto-correlation

% nc = n1_;
% 
% zunit1 = zeros(1,length(round(nc))+5);
% zunit1(round(nc)) = 1;
% acr = xcorr(zunit1,0.2*sr);
% acr(length(acr)/2+0.5) = [];
% acr = reshape(acr,length(acr)/100,100);
% sacr = sum(acr);
% figure;
% bar(sacr)

%% Place map

fminx = min(x_pos(x_pos>0));
fmaxx = max(x_pos);
fminy = min(y_pos(y_pos>0));
fmaxy = max(y_pos);

xedge = fminx-0.0001:8:fmaxx+8;      % spatial bins
yedge = fminy-0.0001:8:fmaxy+8;
xbins = length(xedge) - 1;
ybins = length(yedge) - 1;

% figure;plot(x_pos_ts,x_pos);figure;plot(x_pos_ts,y_pos)
dx = diff(x_pos);       % leave out the outliers
% fdx = find(abs(dx)>8&abs([dx(2:end) 0])>8&abs(dx-[dx(2:end) 0])>16&(x_pos(2:end)|(y_pos(2:end))));
fdx = find(abs(dx)>8&abs([dx(2:end) 0])>8&abs(dx-[dx(2:end) 0])>16);
x_pos(fdx+1) = [];
y_pos(fdx+1) = [];
x_pos_ts(fdx+1) = [];
dy = diff(y_pos);
fdy = find(abs(dy)>8&abs([dy(2:end) 0])>8&abs(dy-[dy(2:end) 0])>16);
x_pos(fdy+1) = [];
y_pos(fdy+1) = [];
x_pos_ts(fdy+1) = [];
% figure;plot(x_pos_ts,x_pos);figure;plot(x_pos_ts,y_pos)

inxs = x_pos > 0 & y_pos > 0;   % leave out (0;0) points
x_pos2 = x_pos(inxs);
y_pos2 = y_pos(inxs);
pos_ts = x_pos_ts(inxs);

ncu = n4;
len = length(ncu);
ap_xpos = zeros(1,len);
ap_ypos = zeros(1,len);
for k = 1:len    % interpolate action potential position
    tap = ncu(k);
    ind = find(pos_ts<tap,1,'last');
    tlow = pos_ts(ind);
    xlow = x_pos2(ind);
    ylow = y_pos2(ind);
    thigh = pos_ts(ind+1);
    xhigh = x_pos2(ind+1);
    yhigh = y_pos2(ind+1);
    mpl = (tap - tlow) / (thigh - tlow);
    ap_xpos(k) = xlow + (xhigh - xlow) * mpl;
    ap_ypos(k) = ylow + (yhigh - ylow) * mpl;
end

hst = hist3([ap_xpos' ap_ypos'],'Edges',{xedge yedge})';     % spike number in space
hst = hst(1:end-1,1:end-1);     % last bin corresponds to the upper edge value
% figure;imagesc(hst)

thst = zeros(size(hst));        % time spent in each bin
tvhst = zeros(size(hst));       % number of visits
for xb = 1:xbins
    for yb = 1:ybins
        xedlow = xedge(xb);
        xedhigh = xedge(xb+1);
        xlin = valuecrossing(pos_ts,x_pos2,xedlow,'up');
        xhin = valuecrossing(pos_ts,x_pos2,xedhigh,'down');
        xlout = valuecrossing(pos_ts,x_pos2,xedlow,'down');
        xhout = valuecrossing(pos_ts,x_pos2,xedhigh,'up');
        xin = union(xlin,xhin);
        xout = union(xlout,xhout);
        if xout(1) < xin(1)
            xin = [pos_ts(1) xin];
        end
        if xin(end) > xout(end)
            xout = [xout pos_ts(end)];
        end
        if ~isequal(length(xin),length(xout))
            error('Technical error 119.')
        end
                
        yedlow = yedge(yb);
        yedhigh = yedge(yb+1);
        ylin = valuecrossing(pos_ts,y_pos2,yedlow,'up');
        yhin = valuecrossing(pos_ts,y_pos2,yedhigh,'down');
        ylout = valuecrossing(pos_ts,y_pos2,yedlow,'down');
        yhout = valuecrossing(pos_ts,y_pos2,yedhigh,'up');
        yin = union(ylin,yhin);
        yout = union(ylout,yhout);
        if yout(1) < yin(1)
            yin = [pos_ts(1) yin];
        end
        if yin(end) > yout(end)
            yout = [yout pos_ts(end)];
        end
        if ~isequal(length(yin),length(yout))
            error('Technical error 137.')
        end
        
        for k1 = 1:length(xin)
            cxin = xin(k1);
            cxout = xout(find(xout>cxin,1,'first'));
            logix = yin<cxout&yout>cxin;
            cyins = yin(logix);
            cyouts = yout(logix);
            for k2 = 1:length(cyins)
                cyin = cyins(k2);
                cyout = cyouts(k2);
                thst(yb,xb) = thst(yb,xb) + (min(cyout,cxout) - max(cyin,cxin));
                tvhst(yb,xb) = tvhst(yb,xb) + 1;
            end
        end
    end
end
% figure;imagesc(thst);
rhst=hst./thst;
figure;imagesc(rhst)

rhst2 = rhst .* double(thst>0.25) .* double(tvhst>3);       % minimum 3 visits, 0.25 s spent in bin
% figure;imagesc(rhst2)
irhst=interp2(rhst2,5);
figure;pcolor(irhst);shading flat

%% Call czplace.m

for k = 1:9
    str=['[irhst' num2str(k) ',rhst_' num2str(k) ']=czplace2(n' num2str(k) ',x_pos,y_pos,x_pos_ts);'];
    eval(str)
end

for k = 1:9     % plot
    str=['subplot(33' num2str(k) ');pcolor(irhst' num2str(k) ');shading flat'];
    eval(str)
    colorbar
end

%% Place correlation (linear)

rhst1n=rhst_1(~isnan(rhst_1));
rhst2n=rhst_2(~isnan(rhst_2));
rhst3n=rhst_3(~isnan(rhst_3));
rhst4n=rhst_4(~isnan(rhst_4));
rhst5n=rhst_5(~isnan(rhst_5));
rhst6n=rhst_6(~isnan(rhst_6));
rhst7n=rhst_7(~isnan(rhst_7));
rhst8n=rhst_8(~isnan(rhst_8));
rhst9n=rhst_9(~isnan(rhst_9));
R = eye(5);

% 3 - 5
r1 = rhst3n;
r2 = rhst5n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(r1,r2);
R(1,2) = corrmtx(2);
R(2,1) = corrmtx(2);

% 3 - 7
r1 = rhst3n;
r2 = rhst7n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(1,3) = corrmtx(2);
R(3,1) = corrmtx(2);

% 3 - 8
r1 = rhst3n;
r2 = rhst8n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(1,4) = corrmtx(2);
R(4,1) = corrmtx(2);

% 3 - 9
r1 = rhst3n;
r2 = rhst9n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(1,5) = corrmtx(2);
R(5,1) = corrmtx(2);

% 5 - 7
r1 = rhst5n;
r2 = rhst7n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(2,3) = corrmtx(2);
R(3,2) = corrmtx(2);

% 5 - 8
r1 = rhst5n;
r2 = rhst8n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(2,4) = corrmtx(2);
R(4,2) = corrmtx(2);

% 5 - 9
r1 = rhst5n;
r2 = rhst9n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(2,5) = corrmtx(2);
R(5,2) = corrmtx(2);

% 7 - 8
r1 = rhst7n;
r2 = rhst8n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(3,4) = corrmtx(2);
R(4,3) = corrmtx(2);

% 7 - 9
r1 = rhst7n;
r2 = rhst9n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(3,5) = corrmtx(2);
R(5,3) = corrmtx(2);

% 8 - 9
r1 = rhst8n;
r2 = rhst9n;
rr1 = r1(r1>0|r2>0);
rr2 = r2(r1>0|r2>0);
corrmtx = corrcoef(rr1,rr2);
R(4,5) = corrmtx(2);
R(5,4) = corrmtx(2);
disp(R)

%% Place correlation - zero pixels included ("Place Field Similarity")

rhs1n=rhs_1(~isnan(rhs_1));
rhs2n=rhs_2(~isnan(rhs_2));
rhs3n=rhs_3(~isnan(rhs_3));
rhs4n=rhs_4(~isnan(rhs_4));
rhs5n=rhs_5(~isnan(rhs_5));
rhs6n=rhs_6(~isnan(rhs_6));
rhs7n=rhs_7(~isnan(rhs_7));
rhs8n=rhs_8(~isnan(rhs_8));
rhs9n=rhs_9(~isnan(rhs_9));
R = eye(5);

% 3 - 5
r1 = rhs3n;
r2 = rhs5n;
corrmtx = corrcoef(r1,r2);
R(1,2) = corrmtx(2);
R(2,1) = corrmtx(2);

% 3 - 7
r1 = rhs3n;
r2 = rhs7n;
corrmtx = corrcoef(r1,r2);
R(1,3) = corrmtx(2);
R(3,1) = corrmtx(2);

% 3 - 8
r1 = rhs3n;
r2 = rhs8n;
corrmtx = corrcoef(r1,r2);
R(1,4) = corrmtx(2);
R(4,1) = corrmtx(2);

% 3 - 9
r1 = rhs3n;
r2 = rhs9n;
corrmtx = corrcoef(r1,r2);
R(1,5) = corrmtx(2);
R(5,1) = corrmtx(2);

% 5 - 7
r1 = rhs5n;
r2 = rhs7n;
corrmtx = corrcoef(r1,r2);
R(2,3) = corrmtx(2);
R(3,2) = corrmtx(2);

% 5 - 8
r1 = rhs5n;
r2 = rhs8n;
corrmtx = corrcoef(r1,r2);
R(2,4) = corrmtx(2);
R(4,2) = corrmtx(2);

% 5 - 9
r1 = rhs5n;
r2 = rhs9n;
corrmtx = corrcoef(r1,r2);
R(2,5) = corrmtx(2);
R(5,2) = corrmtx(2);

% 7 - 8
r1 = rhs7n;
r2 = rhs8n;
corrmtx = corrcoef(r1,r2);
R(3,4) = corrmtx(2);
R(4,3) = corrmtx(2);

% 7 - 9
r1 = rhs7n;
r2 = rhs9n;
corrmtx = corrcoef(r1,r2);
R(3,5) = corrmtx(2);
R(5,3) = corrmtx(2);

% 8 - 9
r1 = rhs8n;
r2 = rhs9n;
corrmtx = corrcoef(r1,r2);
R(4,5) = corrmtx(2);
R(5,4) = corrmtx(2);
disp(R)

%% Place map normalization

rh_1 = rhs_1 / b_max_nonnan(rhs_1);
rh_2 = rhs_2 / b_max_nonnan(rhs_2);
rh_3 = rhs_3 / b_max_nonnan(rhs_3);
rh_4 = rhs_4 / b_max_nonnan(rhs_4);
rh_5 = rhs_5 / b_max_nonnan(rhs_5);
rh_6 = rhs_6 / b_max_nonnan(rhs_6);
rh_7 = rhs_7 / b_max_nonnan(rhs_7);
rh_8 = rhs_8 / b_max_nonnan(rhs_8);
rh_9 = rhs_9 / b_max_nonnan(rhs_9);
rh1n = rh_1(~isnan(rh_1));
rh2n = rh_2(~isnan(rh_2));
rh3n = rh_3(~isnan(rh_3));
rh4n = rh_4(~isnan(rh_4));
rh5n = rh_5(~isnan(rh_5));
rh6n = rh_6(~isnan(rh_6));
rh7n = rh_7(~isnan(rh_7));
rh8n = rh_8(~isnan(rh_8));
rh9n = rh_9(~isnan(rh_9));

%% Cross-correlation

nc1 = n3_;
nc2 = n5_;
spindex = 20;

zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
ccr = xcorr(zunit2,zunit1,0.05*sr);     % 1->2
% ccr(length(acr)/2+0.5) = [];
% ccr = reshape(acr,length(acr)/100,100);
% sacr = sum(acr);
ccr = ccr / length(nc1);
% figure;
% bar(sacr)
% bar(ccr)

% ISI shuffle
leneeg = length(zunit2);
isi = diff(nc2);

lit = length(isi);
rp = randperm(lit);
% while any(rp==1:lit)
%     rp = randperm(lit);
% end
psi1 = [];
for it = 1:lit
    psi1 = [psi1 isi(rp(it))];
end
psvd = [nc2(1) nc2(1)+cumsum(psi1)];
pzunit = zeros(1,length(round(psvd))+5);
pzunit(round(psvd)) = 1;

pccr = xcorr(pzunit,zunit1,0.05*sr);
pccr = pccr / length(nc1);
% figure;
% bar(pccr)

nccr = ccr - pccr;
% figure;
subplot(5,5,spindex)
bar(linspace(-50,50,length(nccr)),nccr)
set(gca,'XLim',[-50 50])

mnc = mean(nccr);
sdc = std(nccr);
thr = mnc + 3 * sdc;
line([-50 50],[thr thr],'Color','red')
line([-50 50],[-thr -thr],'Color','red')

%% Correlation map

rm1 = rhst_3;
rm2 = rhst_5;

s1 = size(rm1,1);
s2 = size(rm1,2);
corrmap = zeros(size(rm1));
for x = 1:s1
    for y = 1:s2
        if ~isnan(rm1(x,y)) && ~isnan(rm2(x,y))
            xind1 = max(x-5,1);
            xind2 = min(x+5,s1);
            yind1 = max(y-5,1);
            yind2 = min(y+5,s2);
            c1 = rm1(xind1:xind2,yind1:yind2);
            c2 = rm2(xind1:xind2,yind1:yind2);
            c11 = c1(~isnan(c1));
            c21 = c2(~isnan(c2));
            C = corrcoef(c11,c21);
            corrmap(x,y) = C(2);
        else
            corrmap(x,y) = NaN;
        end
    end
end
            
figure
pcolor(corrmap)

%% Cross-map

rm1 = rhst_1;
rm2 = rhst_4;

srm1 = (rm1 - b_mean_nonnan(rm1(:))) / b_std_nonnan(rm1(:));
srm2 = (rm2 - b_mean_nonnan(rm2(:))) / b_std_nonnan(rm2(:));

crossmap = srm1 .* srm2;
% figure
% pcolor(crossmap)

icmp = interp2(crossmap,5);
figure
pcolor(icmp)
shading flat

%% Place field with NaNs

for k = 1:9
    str=['rm=rhst_' num2str(k) ';'];
    eval(str)
    rm2 = rm .* zero2nan(double(thst>0.25)) .* zero2nan(double(tvhst>3));
    str=['rhs_' num2str(k) '=rm2;'];
    eval(str)
end

%% KL-distance

% Smoothing
gk = gausskernel([3 3]);
srhs_1 = smooth2_nonnan(rhs_1,gk);
srhs_2 = smooth2_nonnan(rhs_2,gk);
srhs_3 = smooth2_nonnan(rhs_3,gk);
srhs_4 = smooth2_nonnan(rhs_4,gk);
srhs_5 = smooth2_nonnan(rhs_5,gk);
srhs_6 = smooth2_nonnan(rhs_6,gk);
srhs_7 = smooth2_nonnan(rhs_7,gk);
srhs_8 = smooth2_nonnan(rhs_8,gk);
srhs_9 = smooth2_nonnan(rhs_9,gk);

% Normalization (spatial distribution)
rhh_1 = srhs_1 / b_sum_nonnan(srhs_1(:));
rhh_2 = srhs_2 / b_sum_nonnan(srhs_2(:));
rhh_3 = srhs_3 / b_sum_nonnan(srhs_3(:));
rhh_4 = srhs_4 / b_sum_nonnan(srhs_4(:));
rhh_5 = srhs_5 / b_sum_nonnan(srhs_5(:));
rhh_6 = srhs_6 / b_sum_nonnan(srhs_6(:));
rhh_7 = srhs_7 / b_sum_nonnan(srhs_7(:));
rhh_8 = srhs_8 / b_sum_nonnan(srhs_8(:));
rhh_9 = srhs_9 / b_sum_nonnan(srhs_9(:));

rhh1n = rhh_1(~isnan(rhh_1));
rhh2n = rhh_2(~isnan(rhh_2));
rhh3n = rhh_3(~isnan(rhh_3));
rhh4n = rhh_4(~isnan(rhh_4));
rhh5n = rhh_5(~isnan(rhh_5));
rhh6n = rhh_6(~isnan(rhh_6));
rhh7n = rhh_7(~isnan(rhh_7));
rhh8n = rhh_8(~isnan(rhh_8));
rhh9n = rhh_9(~isnan(rhh_9));
D = zeros(5,5);

% 3 - 5
P = rhh3n;
Q = rhh5n;
D(1,2) = KLdist(P,Q);
D(2,1) = KLdist(Q,P);

% 3 - 7
P = rhh3n;
Q = rhh7n;
D(1,3) = KLdist(P,Q);
D(3,1) = KLdist(Q,P);

% 3 - 8
P = rhh3n;
Q = rhh8n;
D(1,4) = KLdist(P,Q);
D(4,1) = KLdist(Q,P);

% 3 - 9
P = rhh3n;
Q = rhh9n;
D(1,5) = KLdist(P,Q);
D(5,1) = KLdist(Q,P);

% 5 - 7
P = rhh5n;
Q = rhh7n;
D(2,3) = KLdist(P,Q);
D(3,2) = KLdist(Q,P);

% 5 - 8
P = rhh5n;
Q = rhh8n;
D(2,4) = KLdist(P,Q);
D(4,2) = KLdist(Q,P);

% 5 - 9
P = rhh5n;
Q = rhh9n;
D(2,5) = KLdist(P,Q);
D(5,2) = KLdist(Q,P);

% 7 - 8
P = rhh7n;
Q = rhh8n;
D(3,4) = KLdist(P,Q);
D(4,3) = KLdist(Q,P);

% 7 - 9
P = rhh7n;
Q = rhh9n;
D(3,5) = KLdist(P,Q);
D(5,3) = KLdist(Q,P);

% 8 - 9
P = rhh8n;
Q = rhh9n;
D(4,5) = KLdist(P,Q);
D(5,4) = KLdist(Q,P);
disp(D)

%% Spatial information: KL-distance from the uniform distribution

v = 1 / length(rhh1n);
unidist = ones(length(rhh1n),1) * v;    % uniform distribution
for k = 1:9
    str = ['KLdist(rhh' num2str(k) 'n,unidist)'];
    eval(str)
end


%% JS-divergence

D = zeros(5,5);

% 3 - 5
P = rhh3n;
Q = rhh5n;
D(1,2) = JSdiv(P,Q);
D(2,1) = JSdiv(Q,P);

% 3 - 7
P = rhh3n;
Q = rhh7n;
D(1,3) = JSdiv(P,Q);
D(3,1) = JSdiv(Q,P);

% 3 - 8
P = rhh3n;
Q = rhh8n;
D(1,4) = KLdist(P,Q);
D(4,1) = KLdist(Q,P);

% 3 - 9
P = rhh3n;
Q = rhh9n;
D(1,5) = JSdiv(P,Q);
D(5,1) = JSdiv(Q,P);

% 5 - 7
P = rhh5n;
Q = rhh7n;
D(2,3) = JSdiv(P,Q);
D(3,2) = JSdiv(Q,P);

% 5 - 8
P = rhh5n;
Q = rhh8n;
D(2,4) = JSdiv(P,Q);
D(4,2) = JSdiv(Q,P);

% 5 - 9
P = rhh5n;
Q = rhh9n;
D(2,5) = JSdiv(P,Q);
D(5,2) = JSdiv(Q,P);

% 7 - 8
P = rhh7n;
Q = rhh8n;
D(3,4) = JSdiv(P,Q);
D(4,3) = JSdiv(Q,P);

% 7 - 9
P = rhh7n;
Q = rhh9n;
D(3,5) = JSdiv(P,Q);
D(5,3) = JSdiv(Q,P);

% 8 - 9
P = rhh8n;
Q = rhh9n;
D(4,5) = JSdiv(P,Q);
D(5,4) = JSdiv(Q,P);
disp(D)

%% f-divergence with abs

D = zeros(5,5);

% 3 - 5
P = rhh3n;
Q = rhh5n;
D(1,2) = fdiv(P,Q,'abs');
D(2,1) = fdiv(Q,P,'abs');

% 3 - 7
P = rhh3n;
Q = rhh7n;
D(1,3) = fdiv(P,Q,'abs');
D(3,1) = fdiv(Q,P,'abs');

% 3 - 8
P = rhh3n;
Q = rhh8n;
D(1,4) = fdiv(P,Q,'abs');
D(4,1) = fdiv(Q,P,'abs');

% 3 - 9
P = rhh3n;
Q = rhh9n;
D(1,5) = fdiv(P,Q,'abs');
D(5,1) = fdiv(Q,P,'abs');

% 5 - 7
P = rhh5n;
Q = rhh7n;
D(2,3) = fdiv(P,Q,'abs');
D(3,2) = fdiv(Q,P,'abs');

% 5 - 8
P = rhh5n;
Q = rhh8n;
D(2,4) = fdiv(P,Q,'abs');
D(4,2) = fdiv(Q,P,'abs');

% 5 - 9
P = rhh5n;
Q = rhh9n;
D(2,5) = fdiv(P,Q,'abs');
D(5,2) = fdiv(Q,P,'abs');

% 7 - 8
P = rhh7n;
Q = rhh8n;
D(3,4) = fdiv(P,Q,'abs');
D(4,3) = fdiv(Q,P,'abs');

% 7 - 9
P = rhh7n;
Q = rhh9n;
D(3,5) = fdiv(P,Q,'abs');
D(5,3) = fdiv(Q,P,'abs');

% 8 - 9
P = rhh8n;
Q = rhh9n;
D(4,5) = fdiv(P,Q,'abs');
D(5,4) = fdiv(Q,P,'abs');
disp(D)

%% Bhattacharyya-distance

D = zeros(5,5);

% 3 - 5
P = rhh3n;
Q = rhh5n;
D(1,2) = Bdist(P,Q);
D(2,1) = Bdist(Q,P);

% 3 - 7
P = rhh3n;
Q = rhh7n;
D(1,3) = Bdist(P,Q);
D(3,1) = Bdist(Q,P);

% 3 - 8
P = rhh3n;
Q = rhh8n;
D(1,4) = Bdist(P,Q);
D(4,1) = Bdist(Q,P);

% 3 - 9
P = rhh3n;
Q = rhh9n;
D(1,5) = Bdist(P,Q);
D(5,1) = Bdist(Q,P);

% 5 - 7
P = rhh5n;
Q = rhh7n;
D(2,3) = Bdist(P,Q);
D(3,2) = Bdist(Q,P);

% 5 - 8
P = rhh5n;
Q = rhh8n;
D(2,4) = Bdist(P,Q);
D(4,2) = Bdist(Q,P);

% 5 - 9
P = rhh5n;
Q = rhh9n;
D(2,5) = Bdist(P,Q);
D(5,2) = Bdist(Q,P);

% 7 - 8
P = rhh7n;
Q = rhh8n;
D(3,4) = Bdist(P,Q);
D(4,3) = Bdist(Q,P);

% 7 - 9
P = rhh7n;
Q = rhh9n;
D(3,5) = Bdist(P,Q);
D(5,3) = Bdist(Q,P);

% 8 - 9
P = rhh8n;
Q = rhh9n;
D(4,5) = Bdist(P,Q);
D(5,4) = Bdist(Q,P);
disp(D)

%% Chi-square divergence

D = zeros(5,5);

% 3 - 5
P = rhh3n;
Q = rhh5n;
D(1,2) = chisquarediv(P,Q);
D(2,1) = chisquarediv(Q,P);

% 3 - 7
P = rhh3n;
Q = rhh7n;
D(1,3) = chisquarediv(P,Q);
D(3,1) = chisquarediv(Q,P);

% 3 - 8
P = rhh3n;
Q = rhh8n;
D(1,4) = chisquarediv(P,Q);
D(4,1) = chisquarediv(Q,P);

% 3 - 9
P = rhh3n;
Q = rhh9n;
D(1,5) = chisquarediv(P,Q);
D(5,1) = chisquarediv(Q,P);

% 5 - 7
P = rhh5n;
Q = rhh7n;
D(2,3) = chisquarediv(P,Q);
D(3,2) = chisquarediv(Q,P);

% 5 - 8
P = rhh5n;
Q = rhh8n;
D(2,4) = chisquarediv(P,Q);
D(4,2) = chisquarediv(Q,P);

% 5 - 9
P = rhh5n;
Q = rhh9n;
D(2,5) = chisquarediv(P,Q);
D(5,2) = chisquarediv(Q,P);

% 7 - 8
P = rhh7n;
Q = rhh8n;
D(3,4) = chisquarediv(P,Q);
D(4,3) = chisquarediv(Q,P);

% 7 - 9
P = rhh7n;
Q = rhh9n;
D(3,5) = chisquarediv(P,Q);
D(5,3) = chisquarediv(Q,P);

% 8 - 9
P = rhh8n;
Q = rhh9n;
D(4,5) = chisquarediv(P,Q);
D(5,4) = chisquarediv(Q,P);
disp(D)

%% Product of spatial firing distributions

% 3 - 5
prd = rhh_3 .* rhh_5;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 3 - 7
prd = rhh_3 .* rhh_7;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 3 - 8
prd = rhh_3 .* rhh_8;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 3 - 9
prd = rhh_3 .* rhh_9;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 3 - 7
prd = rhh_3 .* rhh_7;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 5 - 7
prd = rhh_5 .* rhh_7;
figure
pcolor(prd)
colorbar
shading flat
set(gca,'CLim',[0 3*10^(-5)])

% 5 - 8
prd = rhh_5 .* rhh_8;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 5 - 9
prd = rhh_5 .* rhh_9;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 7 - 8
prd = rhh_7 .* rhh_8;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 7 - 9
prd = rhh_7 .* rhh_9;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

% 8 - 9
prd = rhh_8 .* rhh_9;
figure
pcolor(prd)
shading flat
set(gca,'CLim',[0 3*10^(-5)])
colorbar

%% Firing fields

for k = 1:9
    str = ['s=srhs_' num2str(k) ';'];
    eval(str)
    str = ['ms_' num2str(k) '=s.*(zero2nan(double(s>b_mean_nonnan(s)+0*b_std_nonnan(s))));'];
    eval(str)
end

%% Bhattacharyya-index on firing fields

% Normalization
msn_1 = ms_1 / b_sum_nonnan(ms_1(:));
msn_2 = ms_2 / b_sum_nonnan(ms_2(:));
msn_3 = ms_3 / b_sum_nonnan(ms_3(:));
msn_4 = ms_4 / b_sum_nonnan(ms_4(:));
msn_5 = ms_5 / b_sum_nonnan(ms_5(:));
msn_6 = ms_6 / b_sum_nonnan(ms_6(:));
msn_7 = ms_7 / b_sum_nonnan(ms_7(:));
msn_8 = ms_8 / b_sum_nonnan(ms_8(:));
msn_9 = ms_9 / b_sum_nonnan(ms_9(:));

msn1n = msn_1(~isnan(msn_1));
msn2n = msn_2(~isnan(msn_2));
msn3n = msn_3(~isnan(msn_3));
msn4n = msn_4(~isnan(msn_4));
msn5n = msn_5(~isnan(msn_5));
msn6n = msn_6(~isnan(msn_6));
msn7n = msn_7(~isnan(msn_7));
msn8n = msn_8(~isnan(msn_8));
msn9n = msn_9(~isnan(msn_9));

D = zeros(5,5);

% 3 - 5
P = msn_3;
Q = msn_5;
D(1,2) = 1 - Bdist_nonnan(P,Q);
D(2,1) = 1 - Bdist_nonnan(Q,P);

% 3 - 7
P = msn_3;
Q = msn_7;
D(1,3) = 1 - Bdist_nonnan(P,Q);
D(3,1) = 1 - Bdist_nonnan(Q,P);

% 3 - 8
P = msn_3;
Q = msn_8;
D(1,4) = 1 - Bdist_nonnan(P,Q);
D(4,1) = 1 - Bdist_nonnan(Q,P);

% 3 - 9
P = msn_3;
Q = msn_9;
D(1,5) = 1 - Bdist_nonnan(P,Q);
D(5,1) = 1 - Bdist_nonnan(Q,P);

% 5 - 7
P = msn_5;
Q = msn_7;
D(2,3) = 1 - Bdist_nonnan(P,Q);
D(3,2) = 1 - Bdist_nonnan(Q,P);

% 5 - 8
P = msn_5;
Q = msn_8;
D(2,4) = 1 - Bdist_nonnan(P,Q);
D(4,2) = 1 - Bdist_nonnan(Q,P);

% 5 - 9
P = msn_5;
Q = msn_9;
D(2,5) = 1 - Bdist_nonnan(P,Q);
D(5,2) = 1 - Bdist_nonnan(Q,P);

% 7 - 8
P = msn_7;
Q = msn_8;
D(3,4) = 1 - Bdist_nonnan(P,Q);
D(4,3) = 1 - Bdist_nonnan(Q,P);

% 7 - 9
P = msn_7;
Q = msn_9;
D(3,5) = 1 - Bdist_nonnan(P,Q);
D(5,3) = 1 - Bdist_nonnan(Q,P);

% 8 - 9
P = msn_8;
Q = msn_9;
D(4,5) = 1 - Bdist_nonnan(P,Q);
D(5,4) = 1 - Bdist_nonnan(Q,P);
disp(D)

%% Complementarity

C = zeros(5,5);

% 3 - 5
P = msn_3;
Q = msn_5;
C(1,2) = czcmpl(P,Q,510);
C(2,1) = czcmpl(Q,P,510);

% 3 - 7
P = msn_3;
Q = msn_7;
C(1,3) = czcmpl(P,Q,510);
C(3,1) = czcmpl(Q,P,510);

% 3 - 8
P = msn_3;
Q = msn_8;
C(1,4) = czcmpl(P,Q,510);
C(4,1) = czcmpl(Q,P,510);

% 3 - 9
P = msn_3;
Q = msn_9;
C(1,5) = czcmpl(P,Q,510);
C(5,1) = czcmpl(Q,P,510);

% 5 - 7
P = msn_5;
Q = msn_7;
C(2,3) = czcmpl(P,Q,510);
C(3,2) = czcmpl(Q,P,510);

% 5 - 8
P = msn_5;
Q = msn_8;
C(2,4) = czcmpl(P,Q,510);
C(4,2) = czcmpl(Q,P,510);

% 5 - 9
P = msn_5;
Q = msn_9;
C(2,5) = czcmpl(P,Q,510);
C(5,2) = czcmpl(Q,P,510);

% 7 - 8
P = msn_7;
Q = msn_8;
C(3,4) = czcmpl(P,Q,510);
C(4,3) = czcmpl(Q,P,510);

% 7 - 9
P = msn_7;
Q = msn_9;
C(3,5) = czcmpl(P,Q,510);
C(5,3) = czcmpl(Q,P,510);

% 8 - 9
P = msn_8;
Q = msn_9;
C(4,5) = czcmpl(P,Q,510);
C(5,4) = czcmpl(Q,P,510);
disp(C)

% Visualization
% p=nan2zero(P);
% q=nan2zero(Q);
% figure;pcolor(p+q)

%% Place field size

pfsize = zeros(1,9);
for k = 1:9
    str = ['s=rhst_' num2str(k) ';'];
    eval(str)
    str = ['ms_' num2str(k) '=s.*(zero2nan(double(s>b_mean_nonnan(s)+2*b_std_nonnan(s))));'];
    eval(str)
    str = ['mss_' num2str(k) '=s.*(zero2nan(double(s>b_mean_nonnan(s))));'];
    eval(str)
    str = ['pfsize(k)=czpfsize(rhst_' num2str(k) ');'];
    eval(str)
end
