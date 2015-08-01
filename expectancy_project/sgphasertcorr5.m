function sgphasertcorr5
%SGPHASERTCORR5   Correlation between delta/theta phase and reaction time.
%   SGPHASERTCORR5 calculates linear-circular correlation between 0 ms
%   theta and reaction time (rt.) restricted to those trials where delta
%   phase was in the optimal range. Sorted phase vs. rt. line plots and
%   phase histograms with mean reaction time in each phase bin are saved
%   for pooled samples.
%
%   See also SGPHASERTCORR2, SGPHASERTCORR3 and SGPHASERTCORR4.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];     % input directory
resdir = [DATAPATH 'Expectancy\phasert\bpF_13_47\'];        % output directory
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'delta1-3Hz_theta4-7Hz_singlephase_bandpassFilter.mat'];
load(ff)
data = singleEEGHILB;    % dimensions: indiv.; cond.; EEG chan.;
                         % delta phase-sorted delta phase/theta phase/
                         % reaction time; trials

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
H1 = figure;
H2 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1,H2)
        
        saveas(H1,['deltaboundedtheta_deg60 ' titlestr ' bpF_13_47.fig'])    % save
        saveas(H2,['hist_deltaboundedtheta_deg60 ' titlestr ' bpF_13_47.fig'])    % save
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1,H2)

rt = cell(1,13);    % reaction time
spvr = zeros(13,100,3);
for k = 1:13
    rt{k} = squeeze(data(k,dim2,dim3,3,:));
    spvr(k,:,:) = squeeze(data(k,dim2,dim3,[1 2 3],:))';
end

% Histogram for delta phase
edges = -pi:2*pi/6:pi;      % 60 degree phase bins
psphase_delta = spvr(:,:,1);
sphase_delta = psphase_delta(:);
psphase_theta = spvr(:,:,2);
sphase_theta = psphase_theta(:);
psrt = spvr(:,:,3);
srt = psrt(:);
psp = [sphase_delta sphase_theta srt];
spvr = sortrows(psp,1);
mn_phase = dhist(spvr,edges);
fn = find (mn_phase==max(mn_phase));        % restrict to optimal delta
lim1 = edges(fn);
lim2 = edges(fn+1);
inx = find(sphase_delta>lim1&sphase_delta<lim2);
psp = psp(inx,2:3);

figure(H1);     % sorted theta phase vs. rt.
sp = sortrows(psp,1);
sphase = sp(:,1);
srt = sp(:,2);
[b,bint,r,rint,stats] = regress(srt,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);   % correlation
R = sqrt(stats(1));
p = stats(3);
subplot(2,1,1)
plot(sphase)
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['linear-circular R: ' num2str(R) '  p: ' num2str(p)],'Color','red')
subplot(2,1,2)
plot(smooth(srt,'circular'))
title(titlestr)

edges = -pi:2*pi/9:pi;      % phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
figure(H2);
spvr = sp;
[mn_phase mn_rt] = dhist(spvr,edges);
subplot(2,1,1)
plot([cnts cnts+2*pi],[mn_phase mn_phase],'r')
subplot(2,1,2)
plot([cnts cnts+2*pi],[mn_rt mn_rt],'r')
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