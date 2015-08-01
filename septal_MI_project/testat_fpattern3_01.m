function testat_fpattern3_01
%TESTAT_FPATTERN3   Calculates statistics for transfer entropy data.
%   TESTAT_FPATTERN3 runs on the output of TRTRENRUN. Edit code to specify
%   exact directories!
%
%   It runs on transfer entropy results, and addresses the question whether
%   transfer ebtropy is different between bursting and non-bursting groups
%   in anatomically defined categories of identified cells. It plots and
%   saves boxplots showing the transfer entropy value distributions and the
%   results of the statistical comparisons by Wilcoxon ranksum test.
%
%   See also TESTAT and TESTAT_BURST.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_cont3\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_cont3\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_cont3\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_cont3\control\'];
inpdir_te = [DATAPATH 'TE\timeresolved_tauNTEmax\'];
resdir = [DATAPATH 'TE\Stat_cont3_tauNTEmax\'];
clusdir_theta = [DATAPATH 'MI\Cluster\Theta\'];
mm = pwd;

% Files
[files_theta, files_short_theta] = genfiles(inpdir_theta);
[files_noth, files_short_noth] = genfiles(inpdir_noth);
files_short = intersect(files_short_theta, files_short_noth);
sf = length(files_short);

% Load table
ff = [DATAPATH '5UEntropy\cells'];   % load identification data
[tbl0 tbl] = xlsread(ff);

% Statistics
DF_theta_burst = [];
DF_theta_nonburst = [];
NTE_thetaeu_burst = [];
NTE_thetaeu_nonburst = [];
NTE_thetaue_burst = [];
NTE_thetaue_nonburst = [];
for o = 1:sf
    cd(inpdir_theta)    % load MI values
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_real_theta = aIxy;
    
    cd(inpdir_noth)
    inx = find(strcmp(files_short{o},files_short_noth));
    fnoth = files_noth{inx};
    load(fnoth)
    MI_real_noth = aIxy;
    
    cd(inpdir_theta_control)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_control_theta = aIxy;
    
    cd(inpdir_noth_control)
    inx = find(strcmp(files_short{o},files_short_noth));
    fnoth = files_noth{inx};
    load(fnoth)
    MI_control_noth = aIxy;
    
    cd(clusdir_theta)     % load burst data
    load([ftheta(1:end-15) 'CLUSTER'])
    
    if ~isequal(ftheta(1:6),fnoth(1:6))     % load identification data
        error('Technical error 40')
    else
        fname = ftheta(1:6);
    end
    inx = find(strcmp({tbl{:,1}},fname));
    id1 = tbl{inx,2};
    id2 = tbl{inx,3};
    
    fs = findstr(ftheta,'_');       % load transfer entropy
    ff = [inpdir_te ftheta(1:fs(2)-1) '_TE.mat'];
    load(ff)
    
    dsc = 4;
    fmh = @select;
    if median(MI_real_theta) > median(MI_control_theta)     % check if MI > control
        [Wp_theta,Wh_theta] = b_signrank2(MI_real_theta,MI_control_theta,'alpha',0.01);
    else
        Wh_theta = 0;
    end
    if median(MI_real_noth) > median(MI_control_noth)
        [Wp_noth,Wh_noth] = b_signrank2(MI_real_noth,MI_control_noth,'alpha',0.01);
    else
        Wh_noth = 0;
    end
    if Wh_theta && ClusterNumber == 0
        eval(['DF_theta_nonburst(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
        eval(['NTE_thetaeu_nonburst(end+1:end+dsc) = fmh(NTE_theta_eu,dsc);']);
        eval(['NTE_thetaue_nonburst(end+1:end+dsc) = fmh(NTE_theta_ue,dsc);']);
    elseif Wh_theta && ClusterNumber ~= 0
        eval(['DF_theta_burst(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
        eval(['NTE_thetaeu_burst(end+1:end+dsc) = fmh(NTE_theta_eu,dsc);']);
        eval(['NTE_thetaue_burst(end+1:end+dsc) = fmh(NTE_theta_ue,dsc);']);
    end
end

% Boxplot
cd(resdir)
H = mibox(DF_theta_burst,DF_theta_nonburst,NTE_thetaeu_burst,NTE_thetaeu_nonburst,...
    NTE_thetaue_burst,NTE_thetaue_nonburst,'theta-bursting','non-bursting');
saveas(H,'TEfpattern3_01')
cd(mm)



% -------------------------------------------------------------------------
function [files, files_short] = genfiles(inpdir)

mm = pwd;
cd(inpdir)
files0 = dir(pwd);
files = {};
files_short = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if isequal(ext,'.mat')
        files{end+1} = files0(i).name;
        files_short{end+1} = files0(i).name(1:6);
    end
end
cd(mm)

% -------------------------------------------------------------------------
function H = mibox(m1,m2,p1,p2,q1,q2,l1,l2)

m1 = m1(~isnan(m1)&p1<=1&q1<=1);
m2 = m2(~isnan(m2)&p2<=1&q2<=1);
H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{l1} {l2}],...
    'positions',[0.5 1.1]);
[Wp_PV,Wh_PV] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh_PV
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_PV),'Color',clr,'Horizontalalignment','center')

% -------------------------------------------------------------------------
function X_ds = select(X,dsc)

X_ds = X(1:dsc);