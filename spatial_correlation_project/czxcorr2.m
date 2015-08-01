function [H1 H2] = czxcorr2(ncc1,ncc2)
%CZXCORR2   Crosscorrelation.
%   CZXCORR(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-5 ms time window. Crosscorrelogram is
%   normalized with a shuffled ISI crosscorrelogram to remove random
%   correlations. A significance level of p=0.0013 is indicated on the
%   result. Original and normalized crosscorrelograms are plotted.
%
%   [H1 H2] = CZXCORR2(VD1,VD2) returns the handles of the resulting plots.
%
%   See also XCORR and CZACORR.

% Input argument check
error(nargchk(2,2,nargin))

% Calculate spike times in milliseconds
sr = 2000;  % because it is fed to ccr, it is fine to choose an sr other than 1000
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
nc1(nc1<0.5) = [];
nc2(nc2<0.5) = [];

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
ccr = xcorr(zunit2,zunit1,0.005*sr);     % 1->2; window: -5 ms - 5 ms
ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
H1 = figure;
bar(linspace(-5,5,length(ccr)),ccr)
set(gca,'XLim',[-5 5])

% ISI shuffle
% isi = diff(nc2);
% lit = length(isi);
% rp = randperm(lit);
% psi1 = [];
% for it = 1:lit
%     psi1 = [psi1 isi(rp(it))];
% end
% psvd = [nc2(1) nc2(1)+cumsum(psi1)];
% pzunit = zeros(1,length(round(psvd))+5);
% pzunit(round(psvd)) = 1;
for k = 1:10
    rnd = rand(1) * 100 + 50;
    shf = round(rnd);
    pzunit = [zunit2(shf+1:end) zunit2(1:shf)];
    
% Random crosscorrelogram
    pccr = xcorr(pzunit,zunit1,0.005*sr);
    pccr = pccr / length(nc2);
    str = ['p' num2str(k) '=pccr;'];
    eval(str)
end
pccr = mean([p1;p2;p3;p4;p5;p6;p7;p8;p9;p10]);
% figure;
% bar(pccr)

% Normalized crosscorrelogram
nccr = ccr - pccr;
H2 = figure;
bar(linspace(-5,5,length(nccr)),nccr)
set(gca,'XLim',[-5 5])

mnc = mean(nccr);
sdc = std(nccr);
thr = mnc + 3 * sdc;
nthr = mnc - 3 * sdc;
line([-5 5],[thr thr],'Color','red')      % significance level: p=0.0013
line([-5 5],[nthr nthr],'Color','red')
line([-5 5],[nthr nthr],'Color','red')