function sgphaseampcorr2
%SGPHASEAMPCORR2   Correlation between delta phase and P300 amplitude.
%   SGPHASEAMPCORR2 calculates linear-circular correlation between 0 ms
%   delta phase and P300 amplitude. Amplitude is plotted against sorted
%   deltaphase and regression results are displayed on the plot. Delta
%   phase histograms with mean amplitude in each phase bin are also plotted
%   and saved.
%
%   SGPHASEAMPCORR2 excludes data where the amplitude maximum locaction
%   coincides with one of the limits of the interval within the amplitude
%   maximum was calculated. (In these cases there is no global maximum in
%   the open search interval.)
%
%   See also SGPHASEAMPCORR and SGPHASELATCORR2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\FzCzPz_SE\phase_amp2\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'FiveChannel3_phase_amp_lat.mat'];
load(ff)
data = FiveChannel3_phase_amp_lat;    % dimensions: indiv.; cond.; EEG chan.;
                                      % delta phase/P300 amp./P300 latency; trials
wn = [201 501];     % amplitude search window

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
H1 = figure;
H2 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,wn,H1,H2)
        
        saveas(H1,['sorted_deltaphase_amp ' titlestr])    % save
        saveas(H2,['hist_deltaphase_amp ' titlestr])
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,wn,H1,H2)

deltaphase = cell(1,13);    % extract delta phase and amplitude
lat = cell(1,13);
spvr_delta = zeros(13,100,2);
for k = 1:13
    deltaphase{k} = squeeze(data(k,dim2,dim3,1,:));
    lat{k} = squeeze(data(k,dim2,dim3,2,:));
    deltaphase_vs_lat = squeeze(data(k,dim2,dim3,[1 2],:))';
    spvr_delta(k,:,1:2) = sortrows(deltaphase_vs_lat,1);
end

figure(H1);     % linear-circular regression
spvr = spvr_delta;
psphase = spvr(:,:,1);
sphase = psphase(:);
pslat = spvr(:,:,2);
slat = pslat(:);
psp = [sphase slat];
sp = sortrows(psp,1);
sphase = sp(:,1);
slat = sp(:,2);
inx = find(slat==wn(1)|slat==wn(2));
sphase(inx) = [];
slat(inx) = [];
[b,bint,r,rint,stats] = regress(slat,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);
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
[sslat sds] = smooth(slat,'circular');
lssl = length(sslat);
plot(sslat,1:lssl,'k')
hold on
plot(sslat+sds,1:lssl,'Color',[204 204 204]/255)
plot(sslat-sds,1:lssl,'Color',[204 204 204]/255)
hold off
ylim([0 lssl])
title(titlestr)

edges = -pi:2*pi/18:pi;     % phase histograms
cnts = (edges(1:end-1) + edges(2:end)) / 2;
figure(H2);
spvr = sortrows(psp,1);
inx = find(spvr(:,2)==wn(1)|spvr(:,2)==wn(2));
spvr(inx,:) = [];
[mn_phase mn_lat] = dhist(spvr,edges);
subplot(2,1,1)      % conditional mean amplitude: E(amp|phi1<phase<phi2)
plot([cnts cnts+2*pi],[mn_phase mn_phase],'r')
subplot(2,1,2)
plot([cnts cnts+2*pi],[mn_lat mn_lat],'r')
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
function [mn_phase mn_rt] = dhist(spvr,edges)

n = length(edges);
mn_phase = zeros(1,n-1);
mn_rt = zeros(1,n-1);
for k = 2:n
    inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<edges(k));
    mn_phase(k-1) = length(inx);
    mn_rt(k-1) = mean(spvr(inx,2));
end

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