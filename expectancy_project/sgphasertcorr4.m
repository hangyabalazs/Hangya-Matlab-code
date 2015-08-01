function sgphasertcorr4
%SGPHASERTCORR4   Correlation between delta/theta phase and reaction time.
%   SGPHASERTCORR4 calculates linear-circular correlation between 0 ms
%   delta, theta or delta minus theta phase and reaction time (rt.).
%   Sorted phase vs. rt. line plots are saved for pooled samples.
%
%   See also SGPHASERTCORR2, SGPHASERTCORR3 and SGPHASERTCORR5.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];     % input directory
resdir = [DATAPATH 'Expectancy\phase_rt_ms1200\'];    % output directory
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'FiveChannel3_phase_rt.mat'];
load(ff)
data = FiveChannel3_phase_rt;    % dimensions: indiv.; cond.; EEG chan.;
                                 % delta phase-sorted delta phase/theta phase/
                                 % reaction time; trials

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
H1 = figure;
H2 = figure;
H3 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1,H2,H3)
        
        saveas(H1,['sorted_deltaphase2 ' titlestr '.fig'])    % save
%         saveas(H2,['sorted_thetaphase2 ' titlestr '.fig'])
%         saveas(H3,['sorted_phasediff2 ' titlestr '.fig'])
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1,H2,H3)

thetaphase = cell(1,13);
deltaphase = cell(1,13);
phasediff = cell(1,13);
rt = cell(1,13);
spvr_delta = zeros(13,100,2);
spvr_theta = zeros(13,100,2);
spvr_diff = zeros(13,100,2);
for k = 1:13
    deltaphase{k} = squeeze(data(k,dim2,dim3,1,:));    % theta phase
    thetaphase{k} = squeeze(data(k,dim2,dim3,2,:));    % delta phase
    phasediff{k} = deltaphase{k} - thetaphase{k};      % delta - theta phase
    rt{k} = squeeze(data(k,dim2,dim3,3,:));            % reaction time
            
    deltaphase_vs_rt = squeeze(data(k,dim2,dim3,[1 3],:))';
    spvr_delta(k,:,1:2) = sortrows(deltaphase_vs_rt,1);
    thetaphase_vs_rt = squeeze(data(k,dim2,dim3,[2 3],:))';
    spvr_theta(k,:,1:2) = sortrows(thetaphase_vs_rt,1);
    phasediff_vs_rt = [pmpi(deltaphase_vs_rt(:,1)-thetaphase_vs_rt(:,1))...
        thetaphase_vs_rt(:,2)];
    spvr_diff(k,:,1:2) = sortrows(phasediff_vs_rt,1);
end

figure(H1);     % sorted delta phase vs. rt.
spvr = spvr_delta;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psp = [sphase srt];
sp = sortrows(psp,1);
sphase = sp(:,1);
srt = sp(:,2);
[b,bint,r,rint,stats] = regress(srt,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);      % correlation
R = sqrt(stats(1));
p = stats(3);
subplot(1,2,1)
plot(sphase,1:length(sphase),'k')
x_lim = xlim;
ylim([0 length(sphase)])
y_lim = ylim;
pl = oom(p);
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.1,y_lim(1)+(y_lim(2)-y_lim(1))*0.95,...
    ['R = ' num2str(R) '  p ' pl],'Color','red')
subplot(1,2,2)
[sslat sds] = smooth(srt,'circular');
lssl = length(sslat);
plot(sslat,1:lssl,'k')
hold on
plot(sslat+sds,1:lssl,'Color',[204 204 204]/255)
plot(sslat-sds,1:lssl,'Color',[204 204 204]/255)
hold off
ylim([0 lssl])
title(titlestr)

figure(H2);     % sorted theta phase vs. rt.
spvr = spvr_theta;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psp = [sphase srt];
sp = sortrows(psp,1);
sphase = sp(:,1);
srt = sp(:,2);
[b,bint,r,rint,stats] = regress(srt,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);      % correlation
R = sqrt(stats(1));
p = stats(3);
subplot(2,1,1)
plot(sphase)
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['linear-circular R: ' num2str(R) '  p: ' num2str(p)],'Color','red')
subplot(2,1,2)
[sslat sds] = smooth(srt,'circular');
plot(sslat)
hold on
plot(sslat+sds,'Color',[204 204 204]/255)
plot(sslat-sds,'Color',[204 204 204]/255)
hold off
title(titlestr)

figure(H3);     % sorted phasediff. vs. rt.
spvr = spvr_diff;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psp = [sphase srt];
sp = sortrows(psp,1);
sphase = sp(:,1);
srt = sp(:,2);
[b,bint,r,rint,stats] = regress(srt,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);      % correlation
R = sqrt(stats(1));
p = stats(3);
subplot(2,1,1)
plot(sphase)
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['linear-circular R: ' num2str(R) '  p: ' num2str(p)],'Color','red')
subplot(2,1,2)
[sslat sds] = smooth(srt,'circular');
plot(sslat)
hold on
plot(sslat+sds,'Color',[204 204 204]/255)
plot(sslat-sds,'Color',[204 204 204]/255)
hold off
title(titlestr)

% -------------------------------------------------------------------------
function [X2 S] = smooth(X,str)

n = 101;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
S = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
            S(k) = std(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
                S(k) = std([X(mod2(k-nn,m):m); X(1:k+nn)]) / sqrt(n);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
                S(k) = std([X(k-nn:m); X(1:mod2(k+nn,m))]) / sqrt(n);
            else
                X2(k) = mean(X(k-nn:k+nn));
                S(k) = std(X(k-nn:k+nn)) / sqrt(n);
            end
        end
end

% -------------------------------------------------------------------------
function r = mod2(a,m)

r = mod(a,m);
if r == 0
    r = m;
end

% -------------------------------------------------------------------------
function dt = pmpi(dt)

dt(dt<-pi) = dt(dt<-pi) + 2 * pi;
dt(dt>pi) = dt(dt>pi) - 2 * pi;

% -------------------------------------------------------------------------
function pl = oom(p)

if isequal(p,0)
    pl = '< 10^{-15}';
elseif p < 0.05 && p > 0.01
    pl = '< 0.05';
elseif p < 0.01 && p > 0.001
    pl = '< 0.01';
elseif p < 0.001 && p > 0.0001
    pl = '< 0.001';
elseif p < 0.0001
    ep = (-3:-1:-16);
    ords = 10 .^ ep;
    iop = find(ords>p,1,'last');
    epp = ep(iop);
    pl = ['< 10^{' num2str(epp) '}'];
else
    pl = ['= ' num2str(p)];
end 