function [Stats_eu, Stats_ue, STeu, STue] = emistat_control
%EMISTAT_CONTROL   Calculates statistics for mutual information data.
%   EMISTAT_CONTROL runs on the output of ENTRYRUN3C. Edit code to specify
%   exact directories!
%
%   It runs on entropy results, and addresses the question whether mutual
%   information is different from control in anatomically defined
%   categories of identified cells. It plots and saves boxplots showing the
%   mutual value distributions and the results of the statistical
%   comparisons by Wilcoxon signrank test.
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
MI_theta = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
    'HCNp',[],'HCNp_control',[],'HCNn',[],'HCNn_control',[],...
    'na',[],'na_control',[]);
MI_noth = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
    'HCNp',[],'HCNp_control',[],'HCNn',[],'HCNn_control',[],...
    'na',[],'na_control',[]);
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
    if ~isempty(id1)        % group MI values to anatomical categories
        eval(['MI_theta.' id1 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
        eval(['MI_theta.' id1 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        eval(['MI_noth.' id1 '(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
        eval(['MI_noth.' id1 '_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
    end
    if ~isempty(id2)
        eval(['MI_theta.' id2 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
        eval(['MI_theta.' id2 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        eval(['MI_noth.' id2 '(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
        eval(['MI_noth.' id2 '_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
    end
    if isempty([id1 id2])
        eval(['MI_theta.na(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
        eval(['MI_theta.na_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        eval(['MI_noth.na(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
        eval(['MI_noth.na_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
    end
end

% Boxplot
cd(resdir)
H = mibox(MI_theta.PVp,MI_theta.PVp_control,'PV+','PV+ control');
saveas(H,'MIcontrol_theta_PVp')
H = mibox(MI_theta.HCNp,MI_theta.HCNp_control,'HCN+','HCN+ control');
saveas(H,'MIcontrol_theta_HCNp')
H = mibox(MI_theta.PVn,MI_theta.PVn_control,'PV-','PV- control');
saveas(H,'MIcontrol_theta_PVn')
H = mibox(MI_theta.HCNn,MI_theta.HCNn_control,'HCN-','HCN- control');
saveas(H,'MIcontrol_theta_HCNn')
H = mibox(MI_theta.na,MI_theta.na_control,'na','na control');
saveas(H,'MIcontrol_theta_na')
H = mibox(MI_noth.PVp,MI_noth.PVp_control,'PV+','PV+ control');
saveas(H,'MIcontrol_noth_PVp')
H = mibox(MI_noth.HCNp,MI_noth.HCNp_control,'HCN+','HCN+ control');
saveas(H,'MIcontrol_noth_HCNp')
H = mibox(MI_noth.PVn,MI_noth.PVn_control,'PV-','PV- control');
saveas(H,'MIcontrol_noth_PVn')
H = mibox(MI_noth.HCNn,MI_noth.HCNn_control,'HCN-','HCN- control');
saveas(H,'MIcontrol_noth_HCNn')
H = mibox(MI_noth.na,MI_noth.na_control,'na','na control');
saveas(H,'MIcontrol_noth_na')
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
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{l1} {l2}]);
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