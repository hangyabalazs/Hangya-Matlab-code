function [Stats_eu, Stats_ue, STeu, STue] = emistat_peak
%EMISTAT_PEAK   Calculates statistics for mutual information data.
%   EMISTAT_PEAK runs on the output of ENTRYRUN3C. Edit code to specify
%   exact directories!
%
%   It runs on entropy results, and addresses the question whether mutual
%   information is different among anatomically defined categories of
%   identified cells. It plots and saves a boxplot showing the mutual value
%   distributions and the results of the statistical comparisons by
%   Mann-Whitney U-test.
%
%   EMISTAT_PEAK selects the highest MI values for pooled samples.
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
    if ~isempty(id1)        % group MI values to anatomical categories
        if Wh_theta
            eval(['MI_theta.' id1 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta.' id1 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        end
        if Wh_noth
            eval(['MI_noth.' id1 '(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
            eval(['MI_noth.' id1 '_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
        end
    end
    if ~isempty(id2)
        if Wh_theta
            eval(['MI_theta.' id2 '(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta.' id2 '_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        end
        if Wh_noth
            eval(['MI_noth.' id2 '(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
            eval(['MI_noth.' id2 '_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
        end
    end
    if isempty([id1 id2])
        if Wh_theta
            eval(['MI_theta.na(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
            eval(['MI_theta.na_control(end+1:end+dsc) = fmh(MI_control_theta,dsc);']);
        end
        if Wh_noth
            eval(['MI_noth.na(end+1:end+dsc) = fmh(MI_real_noth,dsc);']);
            eval(['MI_noth.na_control(end+1:end+dsc) = fmh(MI_control_noth,dsc);']);
        end
    end
    
end

% Boxplot
H = mibox(MI_theta);

% Save
cd(resdir)
saveas(H,'MIbox_peak')
save('MIstat_peak','MI_theta','MI_noth')
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
function H = mibox(B)

m1 = B.PVp;
m2 = B.PVn;
m3 = B.HCNp;
m4 = B.HCNn;
m5 = B.na;

H = figure;
boxplot([m1 m2 m3 m4 m5],[zeros(size(m1)) ones(size(m2)) 2*ones(size(m3))...
    3*ones(size(m4)) 4*ones(size(m5))],'labels',[{'PV+'} {'PV-'} {'HCN+'}...
    {'HCN-'} {'na'}]);
[Wp_PV,Wh_PV] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh_PV
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_PV),'Color',clr,'Horizontalalignment','center')
[Wp_HCN,Wh_HCN] = b_ranksum2(m3,m4,'alpha',0.05);
if Wh_HCN
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * 3 / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_HCN),'Color',clr,'Horizontalalignment','center')

% -------------------------------------------------------------------------
function X_ds = select(X,dsc)

sX = sort(X,'descend');
X_ds = sX(1:4);