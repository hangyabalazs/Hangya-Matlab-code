function sgphaseampcorr
%SGPHASEAMPCORR   Correlation between delta phase and P300 amplitude.
%   SGPHASEAMPCORR calculates linear-circular correlation between 0 ms
%   delta phase and P300 amplitude. Amplitude is plotted against sorted
%   deltaphase and regression results are displayed on the plot. Delta
%   phase histograms with mean amplitude in each phase bin are also plotted
%   and saved.
%
%   See also SGPHASEAMPCORR2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\phaseamp3\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'dphase_damp_dlat_new.mat'];
load(ff)
data = dphase_damp_dlat_new;    % dimensions: indiv.; cond.; EEG chan.;
                                % delta phase/P300 amp./P300 latency; trials

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
H1 = figure;
H2 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1,H2)
        
        saveas(H1,['sorted_deltaphase_amp ' titlestr])    % save
        saveas(H2,['hist_deltaphase_amp ' titlestr])
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1,H2)

deltaphase = cell(1,13);    % extract delta phase and amp.
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
[b,bint,r,rint,stats] = regress(slat,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);
R = sqrt(stats(1));
p = stats(3);
subplot(2,1,1)
plot(sphase)
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['linear-circular R: ' num2str(R) '  p: ' num2str(p)],'Color','red')
subplot(2,1,2)
plot(smooth(slat,'circular'))
title(titlestr)

edges = -pi:2*pi/18:pi;     % phase histograms
cnts = (edges(1:end-1) + edges(2:end)) / 2;
figure(H2);
spvr = sortrows(psp,1);
[mn_phase mn_lat] = dhist(spvr,edges);
subplot(2,1,1)
plot([cnts cnts+2*pi],[mn_phase mn_phase],'r')
subplot(2,1,2)      % conditional mean amplitude: E(amp|phi1<phase<phi2)
plot([cnts cnts+2*pi],[mn_lat mn_lat],'r')
title(titlestr)

% -------------------------------------------------------------------------
function X2 = smooth(X,str)

n = 101;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
            else
                X2(k) = mean(X(k-nn:k+nn));
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