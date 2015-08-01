function [Stats_eu, Stats_ue, STeu, STue] = emistat_fpattern2
%EMISTAT_FPATTERN2   Calculates statistics for mutual information data.
%   EMISTAT_FPATTERN2 runs on the output of ENTRYRUN3C. Edit code to specify
%   exact directories!
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
inpdir_theta = [DATAPATH 'MI\Theta\real\'];
inpdir_noth = [DATAPATH 'MI\Noth\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth\control\'];
resdir = [DATAPATH 'MI\Stat\'];
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
MI_theta_burst = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
    'HCNp',[],'HCNp_control',[],'HCNn',[],'HCNn_control',[],...
    'na',[],'na_control',[]);
MI_theta_nonburst = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
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
    
    dsc = 4;
    fmh = @select;
    if ~isempty(id1)        % group MI values to anatomical categories
        if ClusterNumber == 0
            eval(['MI_theta_nonburst.' id1 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta_nonburst.' id1 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        elseif ClusterNumber ~= 0
            eval(['MI_theta_burst.' id1 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta_burst.' id1 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        end
    end
    if ~isempty(id2)
        if ClusterNumber == 0
            eval(['MI_theta_nonburst.' id2 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta_nonburst.' id2 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        elseif ClusterNumber ~= 0
            eval(['MI_theta_burst.' id2 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta_burst.' id2 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        end
    end
    if isempty([id1 id2])
        if ClusterNumber == 0
            eval(['MI_theta_nonburst.na(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta_nonburst.na_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        elseif ClusterNumber ~= 0
            eval(['MI_theta_burst.na(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta_burst.na_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        end
    end
end

% Boxplot
MI_burst = [MI_theta_burst.PVp MI_theta_burst.PVn MI_theta_burst.HCNp ...
    MI_theta_burst.HCNn MI_theta_burst.na];
MI_nonburst = [MI_theta_nonburst.PVp MI_theta_nonburst.PVn MI_theta_nonburst.HCNp ...
    MI_theta_nonburst.HCNn MI_theta_nonburst.na];
cd(resdir)
H = mibox(MI_burst,MI_nonburst,'theta-bursting','non-bursting');
saveas(H,'MIfpattern2')
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