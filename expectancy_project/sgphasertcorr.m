function sgphasertcorr
%SGPHASERTCORR   Correlation between delta/theta phase and reaction time.
%   SGPHASERTCORR calculates linear-circular correlation between 0 ms delta
%   or theta phase and reaction time (rt.). Phase vs. rt. scatter plots and
%   sorted phase vs. rt. line plots are saved for each subject.
%
%   See also SGPHASERTCORR3, SGPHASERTCORR4 and SGPHASERTCORR5.

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
H3 = figure;
H4 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1,H2,H3,H4)
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1,H2,H3,H4)

thetaphase = cell(1,13);    % linear-circular regression
deltaphase = cell(1,13);
rt = cell(1,13);
R_delta = zeros(1,13);
R_theta = zeros(1,13);
p_delta = zeros(1,13);
p_theta = zeros(1,13);
for k = 1:13
    deltaphase{k} = squeeze(data(k,dim2,dim3,1,:));     % delta phase
    thetaphase{k} = squeeze(data(k,dim2,dim3,2,:));     % theta phase
    rt{k} = squeeze(data(k,dim2,dim3,3,:));             % reaction time
    [b,bint,r,rint,stats] = regress(rt{k},...
        [ones(length(deltaphase{k}),1) sin(deltaphase{k}) cos(deltaphase{k})]);
    R_delta(k) = sqrt(stats(1));        % correlation
    p_delta(k) = stats(1);
    [b,bint,r,rint,stats] = regress(rt{k},...
        [ones(length(thetaphase{k}),1) sin(thetaphase{k}) cos(thetaphase{k})]);
    R_theta(k) = sqrt(stats(1));
    p_theta(k) = stats(1);
    
    figure(H1);     % delta phase vs. rt. scatter plot
    plot(deltaphase{k},rt{k},'.')
    ylim([50 700])
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
        ['circular-linear correlation: ' num2str(R_delta(k)) '  p: ' ...
        num2str(p_delta(k))],'Color','red')
    saveas(H1,['deltaphase2 ' titlestr ' ' num2str(k) '.fig'])      % save
    
    figure(H2);     % theta phase vs. rt. scatter plot
    plot(thetaphase{k},rt{k},'.')
    ylim([50 700])
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
        ['circular-linear correlation: ' num2str(R_theta(k)) '  p: ' ...
        num2str(p_theta(k))],'Color','red')
    saveas(H2,['thetaphase2 ' titlestr ' ' num2str(k) '.fig'])      % save
    
    deltaphase_vs_rt = squeeze(data(k,dim2,dim3,[1 3],:))';
    spvr = sortrows(deltaphase_vs_rt,1);
    figure(H3);     % sorted delta phase vs. rt.
    subplot(2,1,1)
    plot(spvr(:,1))
    subplot(2,1,2)
    plot(spvr(:,2))
    saveas(H3,['sorted_deltaphase ' titlestr ' ' num2str(k) '.fig'])    % save
    
    thetaphase_vs_rt = squeeze(data(k,dim2,dim3,[2 3],:))';
    spvr = sortrows(thetaphase_vs_rt,1);
    figure(H4);     % sorted theta phase vs. rt.
    subplot(2,1,1)
    plot(spvr(:,1))
    subplot(2,1,2)
    plot(spvr(:,2))
    saveas(H4,['sorted_thetaphase ' titlestr ' ' num2str(k) '.fig'])    % save
end

close all