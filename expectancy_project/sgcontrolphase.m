function sgcontrolphase
%SGCONTROLPHASE   Comparison of delta phase dispersion in control expectancy data.
%   SGCONTROLPHASE compares kappa in the 20% and in the 80% condition of
%   expectancy control data. Note that the because of the unbalanced pooled
%   sample in SGCONTROLPHASE, SGCONTROLPHASE2 should be used instead!
%
%   See also KAPPACOMPARE2 and SGCONTROLPHASE2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\control\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'kontroll_hilbert_deltaphase.mat'];
load(ff)
ng1 = sh21;     % 20% condition
ng2 = sh41;     % 80% condition

% Circular statistics
[ftm1, mn1, mvl1] = mvlmn(ng1,'rad');      % first trigonometric moment, circ. mean and mean res. length
[kappa_ML1 kappa1] = kappaest(ng1,'rad');  % concentration parameter
[ftm2, mn2, mvl2] = mvlmn(ng2,'rad');
[kappa_ML2 kappa2] = kappaest(ng2,'rad');
[p, H] = kappacompare2(ng1,ng2,'rad','twosided');   % Monte Carlo randomisation test

% Plot
edges2 = -pi:2*pi/9:pi;     % bin limits for phase hist.
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
H = figure;     % phase histograms (blue: 20%, red: 80%)
subplot(1,2,1)
hold on
[nm xout] = histc(ng1,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm nm],'b')
ylim([40 180])
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    ['mean: ' num2str(mn1)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    ['mvl: ' num2str(mvl1)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    ['kappa: ' num2str(kappa1)],'Color','blue')
title('20%')
subplot(1,2,2)
hold on
[nm xout] = histc(ng2,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm nm],'r')
ylim([40 180])
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
title('80%')

% Save
saveas(H,'control_kappacomp.fig')
cd(mm)