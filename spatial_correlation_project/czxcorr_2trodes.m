function czxcorr_2trodes(rhst1,rhst2,vdisc1,vdisc2)
%CZXCORR_2TRODES   Cross-correlation for cells recorded from different trodes.
%   CZXCORR_2TRODES(RH1,RH2,VD_PYR,VD_INT) calculates complementarity index 
%   for rate maps (RH1: prramidal cell, RH2: interneuron), linear
%   correlation and normalized crosscorrelation for a pyramidal cell
%   (VD_PYR) and an interneuron (VD_INT).
%
%   See also CZCMPL2, CZXCORR and CZXCORR2.

% Directories
intcode = 'nr011';
pyrcode = 'nr001';
global DATAPATH
resdir = [DATAPATH 'Czurko\discriminated2\new2\acin08s045l\Recluster_' ...
    intcode '_' pyrcode '\'];
if ~isdir(resdir)
    mkdir(resdir)
end

% Complementarity index, linear correlation (place field similarity)
s_int = rhst2;
ms3_int = s_int .* (zero2nan(double(s_int>b_max_nonnan(s_int)*0.5)));
s_pyr = rhst1;
ms3_pyr = s_pyr .* (zero2nan(double(s_pyr>b_max_nonnan(s_pyr)*0.5)));
c = length(s_int(~isnan(s_int)));
[C Cb Cc] = czcmpl2(ms3_int,ms3_pyr,c);
R = czpfs(s_int,s_pyr);
Rmod = czpfs_mod(s_int,s_pyr);
fn = [resdir 'complementarity_' intcode '_' pyrcode '.mat'];       % result file
save(fn,'C','Cb','Cc','R','Rmod')

% Normalized cross-correlation
vdisc_x = vdisc1;
vdisc_y = vdisc2;
[H1 H2] = lczxcorr(vdisc_x',vdisc_y');   % window: +-50 ms
figure(H1)
A_cross = gca;
ach1 = findobj(allchild(A_cross),'Type','line');    % crosscorrelation
ach2 = findobj(allchild(A_cross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
set(A_cross,'XLim',[-50 50])
xd = get(ach2,'XData');
yd = get(ach2,'YData');
axes(A_cross)
axis off
fns = [resdir 'CROSS_' intcode '_' pyrcode '.fig'];
saveas(gcf,fns)    % save
figure(H2)
A_normcross = gca;
ach1 = findobj(allchild(A_normcross),'Type','line');    % crosscorrelation
ach2 = findobj(allchild(A_normcross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
yd = get(ach2,'YData');
trsr = yd(52) + yd(53) + yd(54);
set(A_normcross,'XLim',[-50 50])
xd = get(ach2,'XData');
yd = get(ach2,'YData');
xdl = get(ach1(3),'XData');
ydl = get(ach1(3),'YData');
axes(A_normcross)
axis off
fns = [resdir 'NCROSS_' intcode '_' pyrcode '.fig'];
saveas(gcf,fns)    % save
fn = [resdir 'trsr_' intcode '_' pyrcode '.mat'];       % result file
save(fn,'trsr')

[H1 H2] = lczxcorr2(vdisc_x',vdisc_y');   % window: +-5 ms
figure(H1)
A_smallcross = gca;
ach1 = findobj(allchild(A_smallcross),'Type','line');    % crosscorrelation
ach2 = findobj(allchild(A_smallcross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
set(A_smallcross,'XLim',[-5 5])
xd = get(ach2,'XData');
yd = get(ach2,'YData');
axes(A_smallcross)
axis off
fns = [resdir 'SCROSS_' intcode '_' pyrcode '.fig'];
saveas(gcf,fns)    % save

% -------------------------------------------------------------------------
function [H1 H2 trsc] = lczxcorr(ncc1,ncc2)
%CZXCORR   Crosscorrelation.
%   CZXCORR(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-50 ms time window. Crosscorrelogram is
%   normalized with a shuffled ISI crosscorrelogram to remove random
%   correlations. A significance level of p=0.0013 is indicated on the
%   result. Original and normalized crosscorrelograms are plotted.
%
%   [H1 H2] = CZXCORR(VD1,VD2) returns the handles of the resulting plots.
%
%   [H1 H2 TRSC] = CZXCORR(VD1,VD2) returns transmission success rate as
%   well. Note, that interneuron has to be the second input argument and 
%   pyramidal cell has to be the first input argument for transmission
%   success rate calculation.
%
%   See also XCORR and CZACORR.

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
ccr = xcorr(zunit2,zunit1,0.05*sr);     % 1->2; window: -50 ms - 50 ms
ccr = ccr / length(nc1);     % norm. with the no. of ref. events to get transmission prob.
H1 = figure;
bar(linspace(-50,50,length(ccr)),ccr)
set(gca,'XLim',[-50 50])

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
    pccr = xcorr(pzunit,zunit1,0.05*sr);
    pccr = pccr / length(nc1);
    str = ['p' num2str(k) '=pccr;'];
    eval(str)
end
pccr = mean([p1;p2;p3;p4;p5;p6;p7;p8;p9;p10]);
% figure;
% bar(pccr)

% Normalized crosscorrelogram
nccr = ccr - pccr;
trsc = max(nccr(45:55));      % transmission success rate
H2 = figure;
bar(linspace(-50,50,length(nccr)),nccr)
set(gca,'XLim',[-50 50])

mnc = mean(nccr);
sdc = std(nccr);
thr = mnc + 3 * sdc;
nthr = mnc - 3 * sdc;
line([-50 50],[thr thr],'Color','red')      % significance level: p=0.0013
line([-50 50],[nthr nthr],'Color','red')
line([-50 50],[nthr nthr],'Color','red')

% -------------------------------------------------------------------------
function [H1 H2] = lczxcorr2(ncc1,ncc2)
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
ccr = ccr / length(nc1);     % norm. with the no. of ref. events to get transmission prob.
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
    pccr = pccr / length(nc1);
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
