function [Stats_eu, Stats_ue, STeu, STue] = emistat_control_hc
%EMISTAT_CONTROL_HC   Calculates statistics for mutual information data.
%   EMISTAT_CONTROL_HC runs on the output of ENTRYRUN3C. Edit code to
%   specify exact directories!
%
%   It runs on entropy results, and addresses the question whether mutual
%   information is different from control. It plots and saves a boxplot
%   showing the mutual value distributions and the results of the
%   statistical comparisons by Wilcoxon signrank test.
%
%   See also ENTRYSTAT6.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_hc\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_hc\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_hc\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_hc\control\'];
resdir = [DATAPATH 'MI\Stat_hc\'];
mm = pwd;

% Files
[files_theta, files_short_theta] = genfiles(inpdir_theta);
[files_noth, files_short_noth] = genfiles(inpdir_noth);
files_short = intersect(files_short_theta, files_short_noth);
sf = length(files_short);

% Statistics
MI_theta = [];
MI_theta_control = [];
MI_noth = [];
MI_noth_control = [];
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
    
    dsc = 4;
    fmh = @select;
    eval(['MI_theta(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
    eval(['MI_theta_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
    eval(['MI_noth(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
    eval(['MI_noth_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
end

% Boxplot
H1 = mibox(MI_theta,MI_theta_control,'HC theta','HC theta control');
H2 = mibox(MI_noth,MI_noth_control,'HC noth','HC noth control');

% Save
cd(resdir)
saveas(H1,'MIbox_theta_control_hc')
saveas(H2,'MIbox_noth_control_hc')
save('MIstat','MI_theta','MI_noth')
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