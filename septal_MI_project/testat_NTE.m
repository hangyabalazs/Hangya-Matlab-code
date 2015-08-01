function testat_NTE
%TESTAT_NTE   Calculates statistics for transfer entropy data.
%   TESTAT_NTE runs on the output of TRTRENRUN. Edit code to specify exact
%   directories!
%
%   It runs on transfer entropy results, and addresses the question whether
%   MS->HC and HC->MS normalized transfer entropy (NTE) values are
%   different during theta and non-theta segments. It plots and saves
%   boxplots showing the NTE value distributions and the results of the
%   statistical comparisons by Wilcoxon ranksum test.
%
%   See also TESTAT, TESTAT_BURST and TESTAT_FPATTERN3.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_cont3\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_cont3\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_cont3\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_cont3\control\'];
inpdir_te = [DATAPATH 'TE\timeresolved_tauNTEmax\'];
resdir = [DATAPATH 'TE\Stat_cont3_tauNTEmax\'];
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
MSHC_theta = [];
HCMS_theta = [];
MSHC_noth = [];
HCMS_noth = [];
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
    if median(MI_real_theta) > median(MI_real_noth)
        Wh = Wh_theta;
    else
        Wh = Wh_noth;
    end
    if Wh
        eval(['HCMS_theta(end+1:end+dsc) = fmh(NTE_theta_eu,dsc);']);
        eval(['MSHC_theta(end+1:end+dsc) = fmh(NTE_theta_ue,dsc);']);
    end
    if Wh
        eval(['HCMS_noth(end+1:end+dsc) = fmh(NTE_noth_eu,dsc);']);
        eval(['MSHC_noth(end+1:end+dsc) = fmh(NTE_noth_ue,dsc);']);
    end
end

% Boxplot
cd(resdir)
H = mibox(HCMS_theta,MSHC_theta,HCMS_noth,MSHC_noth,...
    'HC->MS theta','MS->HC theta','HC->MS noth','MS->HC noth');
saveas(H,'TE_NTE01')
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
function H = mibox(pm1,pm2,pm3,pm4,l1,l2,l3,l4)

m1 = pm1((pm1>0|pm2>0)&(pm1<=1&pm2<=1));
m2 = pm2((pm1>0|pm2>0)&(pm1<=1&pm2<=1));
m3 = pm3((pm3>0|pm4>0)&(pm3<=1&pm4<=1));
m4 = pm4((pm3>0|pm4>0)&(pm3<=1&pm4<=1));
H = figure;
boxplot([m1 m2 m3 m4],[zeros(size(m1)) ones(size(m2)) 2*ones(size(m3)) 3*ones(size(m4))],...
    'labels',[{l1} {l2} {l3} {l4}],'positions',[1 1.8 3.3 4.1]);
[Wp_PV1,Wh_PV1] = b_signrank2(m1,m2,'alpha',0.05);
[Wp_PV2,Wh_PV2] = b_signrank2(m3,m4,'alpha',0.05);
if Wh_PV1
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
y_lim = [0 2];
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_PV1),'Color',clr,'Horizontalalignment','center')
if Wh_PV2
    clr = 'red';
else
    clr = 'black';
end
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * 3 / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_PV2),'Color',clr,'Horizontalalignment','center')

% -------------------------------------------------------------------------
function X_ds = select(X,dsc)

X_ds = X(1:dsc);