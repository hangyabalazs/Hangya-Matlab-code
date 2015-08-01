function sgcontrolphase2
%SGCONTROLPHASE2   Comparison of delta phase dispersion in control expectancy data.
%   SGCONTROLPHASE2 compares kappa in the 20% and in the 80% condition of
%   expectancy control data. Indiviual phase distributions are overlayed on
%   the first plot, whereas a pooled sample is generated to serve a basis
%   for the statistics of the second plot, choosing a random sample of each
%   individual distribution of the same size. The sample size is determined
%   by the smallest individual sample (371 and 89 for the 20% and 80% 
%   condition, respectively). Monte Carlo permutation test is applied to 
%   compare the concentration parameter of the two samples. The test is
%   built on the concentration parameter itself and performs 10000
%   randomisation step to estimate the significance level (See
%   KAPPACOMPARE2 for details.)
%
%   See also KAPPACOMPARE2 and SGCONTROLPHASE2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\expectancy_rawDATA_EXP2nofilt\'];
resdir = [DATAPATH 'Expectancy\filtart_control_filtafterepoch_gamma\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'filtartgamma'];
load(ff)
xz = 'Fz';      % channel
sh211 = eval(['sh21' upper(xz)]);
sh411 = eval(['sh41' upper(xz)]);

% Plot all individuals
edges2 = -pi:2*pi/18:pi;     % bin limits for phase hist.
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
H1 = figure;    % individual delta phase distributions (blue: 20%, red: 80%)
ng1 = [];
ng2 = [];
ng1_indiv = cell(1,11);
ng2_indiv = cell(1,11);
for k = 1:11
    ngk1 = sh211(k).d';     % 20%
    ngk2 = sh411(k).d';     % 80%
    if size(ngk1,2) == 1
        ngk1 = ngk1';
    end        
    if size(ngk2,2) == 1
        ngk2 = ngk2';
    end        
    subplot(1,2,1)
    hold on
    [nm xout] = histc(ngk1,edges2);
    nm = nm(1:end-1);
    plot([cnts2 cnts2+2*pi],[nm/sum(nm) nm/sum(nm)],'b')
    title('20%')
    subplot(1,2,2)
    hold on
    [nm xout] = histc(ngk2,edges2);
    nm = nm(1:end-1);
    plot([cnts2 cnts2+2*pi],[nm/sum(nm) nm/sum(nm)],'r')
    title('80%')
    
    rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
    rp1 = randperm(length(ngk1));       % downsample
    rp2 = randperm(length(ngk2));
    ng1_indiv{k} = ngk1(rp1(1:371));
    ng2_indiv{k} = ngk2(rp2(1:89));
    ng1 = [ng1 ngk1(rp1(1:371))];       % pooled sample
    ng2 = [ng2 ngk2(rp2(1:89))];
end
    

% Circular statistics
[ftm1, mn1, mvl1] = mvlmn(ng1,'rad');      % first trigonometric moment, circ. mean and mean res. length
[kappa_ML1 kappa1] = kappaest(ng1,'rad');  % concentration parameter
[ftm2, mn2, mvl2] = mvlmn(ng2,'rad');
[kappa_ML2 kappa2] = kappaest(ng2,'rad');
[p, H] = kappacompare2(ng1,ng2,'rad','twosided');   % Monte Carlo rand. test
n = length(ng1);
z = n * (mvl1 ^ 2);  % Rayleigh's Z statistic
pR1 = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
    (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
n = length(ng2);
z = n * (mvl2 ^ 2);
pR2 = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
    (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));

% Summary plot
H2 = figure;     % phase histograms (blue: 20%, red: 80%)
subplot(1,2,1)
hold on
[nm xout] = histc(ng1,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm/sum(nm) nm/sum(nm)],'b')
ylim([0.05 0.18])
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['mean: ' num2str(mn1)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    ['mvl: ' num2str(mvl1)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    ['kappa: ' num2str(kappa1)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.6,...
    ['Rayleigh p : ' num2str(pR1)],'Color','blue')
title('20%')
subplot(1,2,2)
hold on
[nm xout] = histc(ng2,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm/sum(nm) nm/sum(nm)],'r')
ylim([0.05 0.18])
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['mean: ' num2str(mn2)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    ['mvl: ' num2str(mvl2)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    ['kappa: ' num2str(kappa2)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.6,...
    ['p (perm. test): ' num2str(p)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.5,...
    ['Rayleigh p : ' num2str(pR2)],'Color','red')
title('80%')

% Plots for publication
H3 = figure;    % rose diagram
[tt,rr] = rose(ones(1,5));
R = polar(tt,rr/100);
hold on
[tout rout] = rose(ng1);
polar(tout,rout/sum(rout))
hold off
delete(R)
title('20%')
H4 = figure;
[tt,rr] = rose(ones(1,5));
R = polar(tt,rr/100);
hold on
[tout rout] = rose(ng2);
polar(tout,rout/sum(rout))
hold off
delete(R)
title('80%')

H5 = figure;    % phase histogram
for k = 1:11
    [nm xout] = histc(ng1_indiv{k},edges2);
    nm = nm(1:end-1);
    nm = nm / sum(nm);
    plot([cnts2 cnts2+2*pi],[nm nm],'Color',[200 200 200]/255,'LineWidth',1)
    hold on
end
[nm xout] = histc(ng1,edges2);
nm = nm(1:end-1);
nm = nm / sum(nm);      % show probabilities
plot([cnts2 cnts2+2*pi],[nm nm],'k','LineWidth',2)
ylim([0 0.16])
xlim([min(cnts2) max(cnts2+2*pi)])
hold off
H6 = figure;
for k = 1:11
    [nm xout] = histc(ng2_indiv{k},edges2);
    nm = nm(1:end-1);
    nm = nm / sum(nm);
    plot([cnts2 cnts2+2*pi],[nm nm],'Color',[200 200 200]/255,'LineWidth',1)
    hold on
end
[nm xout] = histc(ng2,edges2);
nm = nm(1:end-1);
nm = nm / sum(nm);      % show probabilities
plot([cnts2 cnts2+2*pi],[nm nm],'k','LineWidth',2)
ylim([0 0.16])
xlim([min(cnts2) max(cnts2+2*pi)])
hold off

% Save
saveas(H1,['control_alldistrib_' xz '.fig'])
saveas(H2,['control_kappacomp_' xz '.fig'])
saveas(H3,['control_rose_20%_' xz '.fig'])
saveas(H4,['control_rose_80%_' xz '.fig'])
saveas(H5,['control_phasehist_20%_' xz '.fig'])
saveas(H6,['control_phasehist_80%_' xz '.fig'])
cd(mm)