function [Stats_eu, Stats_ue, STeu, STue] = emistat_fpattern_hc
%EMISTAT_FPATTERN_HC   Calculates statistics for mutual information data.
%   EMISTAT_FPATTERN_HC runs on the output of ENTRYRUN3C. Edit code to
%   specify exact directories!
%
%   It runs on entropy results, and addresses the question whether mutual
%   information is different between bursting and non-bursting groups in
%   anatomically defined categories of identified cells. It plots and saves
%   boxplots showing the mutual value distributions and the results of the
%   statistical comparisons by Wilcoxon ranksum test.
%
%   See also EMISTAT and ENTRYSTAT6.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_hc_cont3\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_hc_cont3\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_hc_cont3\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_hc_cont3\control\'];
resdir = [DATAPATH 'MI\Stat_hc_cont3\'];
clusdir_theta = [DATAPATH 'Burst\Cluster\Theta_hc\'];
mm = pwd;

% Files
[files_theta, files_short_theta] = genfiles(inpdir_theta);
[files_noth, files_short_noth] = genfiles(inpdir_noth);
files_short = intersect(files_short_theta, files_short_noth);
sf = length(files_short);

% Statistics
MI_theta_burst = [];
MI_theta_nonburst = [];
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
        eval(['MI_theta_nonburst(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
    elseif Wh_theta && ClusterNumber ~= 0
        eval(['MI_theta_burst(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
    end
end

% Boxplot
cd(resdir)
H = mibox(MI_theta_burst,MI_theta_nonburst,'theta-bursting','non-bursting');
saveas(H,'MIfpattern_hc_cont3')
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

X_ds = X(1:4);