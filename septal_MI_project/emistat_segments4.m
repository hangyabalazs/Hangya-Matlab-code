function [Stats_eu, Stats_ue, STeu, STue] = emistat_segments3
%EMISTAT_SEGMENTS3   Calculates statistics for mutual information data.
%   EMISTAT_SEGMENTS3 runs on the output of ENTRYRUN3C. Edit code to specify
%   exact directories!
%
%   It runs on entropy results, and addresses the question whether mutual
%   information is different between theta and non-theta segments. It plots
%   and saves boxplots showing the mutual value distributions and the
%   results of the statistical comparisons by Wilcoxon ranksum test.
%
%   See also EMISTAT and ENTRYSTAT6.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta\real\'];
inpdir_noth = [DATAPATH 'MI\Noth\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth\control\'];
resdir = [DATAPATH 'MI\Stat\'];
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
MI_theta = [];
MI_noth = [];
counter_theta = 0;
counter_noth = 0;
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
        if Wh_theta
            [Wp2,Wh2] = b_ranksum2(MI_real_theta,MI_real_noth,'alpha',0.01);
            if Wp2
                counter_theta = counter_theta + 1;
            end
        end
    else
        Wh = Wh_noth;
        if Wh_noth
            [Wp2,Wh2] = b_ranksum2(MI_real_noth,MI_real_theta,'alpha',0.01);
            if Wp2
                counter_noth = counter_noth + 1;
            end
        end
    end
    if Wh
        eval(['MI_theta(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
        eval(['MI_noth(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
    end
end

% Boxplot
cd(resdir)
H = mibox(MI_theta,MI_noth,'theta','noth');
saveas(H,'MIsegments4')
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
function H = mibox(m1,m2,l1,l2)

H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{l1} {l2}],...
    'positions',[0.5 1.1]);
[Wp_PV,Wh_PV] = b_signrank2(m1,m2,'alpha',0.05);
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

X_ds = X(1:4);