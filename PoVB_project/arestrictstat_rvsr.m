function arestrictstat_rvsr
%ARESTRICTSTAT_RVSR   Beta MVL vs. slow oscillation MVL.
%   ARESTRICTSTAT_RVSR calculates mean resultant length for beta and slow
%   oscillation and calculates correlation between them. It also correlates
%   beta and slow osclillation mean angle. Edit code to modify input and
%   output directories!
%
%   See also ARESTRICTSTAT_PHASE3_GAMMA.

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_spindle_stand\'];   % beta phase analysis data
inpdir3 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_stand\'];   % slow osc. phase analysis data
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_spindle_stand\'];
mm = pwd;
dr = dir(inpdir);

% Main
H1 = figure;
hold on
H2 = figure;
hold on
clr.Po = 'r';
clr.VB = 'b';
clr.PoVPM = 'k';
clr.LD = 'm';
clr.nRT = 'g';
cum_slow_mvl = [];
cum_beta_mvl = [];
cum_slow_ang = [];
cum_beta_ang = [];
for o = 3:length(dr)
    inpadd = dr(o).name;    % load burst data
    cd(inpdir)
    cd(inpadd)
    cd('bas')
    ddr = dir(pwd);
    fn = ddr(end).name(1:end-4);
    cmps = strread(fn,'%s','delimiter','_');
    fname = [cmps{1} '_' cmps{2}];
    ff = [inpdir2 fn '_PHASE.mat'];
    try
        load(ff)
    catch
        lasterr
        continue
    end
    ff = [tabledir 'tablazat_Balazsnak2'];   % load position data
    [tbl0 tbl] = xlsread(ff);
    inx = find(strcmp({tbl{:,1}},fname));
    if isempty(inx)
        continue
    end
    loc = tbl{inx,3};
    eval(['clrr = clr.' loc ';']);
    
    n_afsp = length(aang_afsp);     % all first spikes
    z = n_afsp * (mvl_afsp ^ 2);  % Rayleigh's Z statistic
    p_afsp = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
        (4 * n_afsp) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n_afsp ^ 2));
    beta_mvl = mvl_afsp;
    beta_ang = circular_mean(aang_afsp,'deg');
    
    ff = [inpdir3 fn '_PHASE.mat'];     % load slow osc. phase data
    load(ff)
    [ftm, ang, mvl] = mvlmn(aang_afsp,'deg');
    slow_mvl = mvl;
    slow_ang = circular_mean(aang_afsp,'deg');
    
    if p_afsp < 0.01
        figure(H1)
        plot(slow_mvl,beta_mvl,[clrr '.'],'MarkerSize',20);
        figure(H2)
        plot(slow_ang,beta_ang,[clrr '.'],'MarkerSize',20);
        cum_slow_mvl(end+1) = slow_mvl;
        cum_slow_ang(end+1) = slow_ang;
        cum_beta_mvl(end+1) = beta_mvl;
        cum_beta_ang(end+1) = beta_ang;
    end
end

% Plot and save
dbclear if error
cd(resdir)
figure(H1)
xlabel('Slow osc. MVL')
ylabel('beta MVL')
[gr icp err] = linefit(cum_slow_mvl,cum_beta_mvl);  % fit line
x = [min(cum_slow_mvl):0.01:max(cum_slow_mvl)];
y = x .* gr + icp;
hold on
plot(x,y)
[b,bint,r,rint,stats] = regress(cum_slow_mvl',[ones(length(cum_beta_mvl),1) cum_beta_mvl']);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(cum_slow_mvl,cum_beta_mvl);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
saveas(H1,'rvsr_new.fig')
saveas(H1,'rvsr_new.eps')

figure(H2)
xlabel('Slow osc. angle')
ylabel('beta angle')
saveas(H2,'avsa_new.fig')
saveas(H2,'avsa_new.eps')

H3 = figure;
cum_beta_ang2 = cum_beta_ang;
cum_beta_ang2(cum_beta_ang2<0) = cum_beta_ang2(cum_beta_ang2<0) + 360;
plot(cum_slow_ang,cum_beta_ang2,'.')
[gr icp err] = linefit(cum_slow_ang,cum_beta_ang2);  % fit line
x = [min(cum_slow_ang):0.01:max(cum_slow_ang)];
y = x .* gr + icp;
hold on
plot(x,y)
[b,bint,r,rint,stats] = regress(cum_slow_ang',[ones(length(cum_beta_ang2),1) cum_beta_ang2']);
R1 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
preR = corrcoef(cum_slow_ang,cum_beta_ang2);
R = preR(2);
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.9+y_lim(1),['R: ' num2str(R)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.8+y_lim(1),['F: ' num2str(F)])
text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.7+y_lim(1),['p (F-test): ' num2str(p)])
saveas(H3,'avsa2_new.fig')
saveas(H3,'avsa2_new.eps')