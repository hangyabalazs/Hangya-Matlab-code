function testat
%TESTAT   Calculates statistics for transfer entropy data.
%   TESTAT runs on the output of TRTRENRUN. Edit code to specify exact
%   directories!
%
%   It runs on transfer entropy results, and addresses the question whether
%   transfer entropy is different among anatomically defined categories of
%   identified cells. It plots and saves a boxplot showing the transfer
%   entropy value distributions and the results of the statistical
%   comparisons by Mann-Whitney U-test.
%
%   See also TESTAT_BURST and TESTAT_FPATTERN3.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_cont3\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_cont3\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_cont3\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_cont3\control\'];
inpdir_te = [DATAPATH 'TE\timeresolved_tauNTEmax_modTP\'];
resdir = [DATAPATH 'TE\Stat_cont3_tauNTEmax_modTP\'];
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
cellcounter = 0;
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
    DIF_theta = NTE_theta_eu - NTE_theta_ue;
    DIF_noth = NTE_noth_eu - NTE_noth_ue;
    DIF_theta = DF_theta_eu;
    DIF_noth = DF_theta_eu;
    
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
    if Wh_theta
        cellcounter = cellcounter + 1;
    end
    if ~isempty(id1)        % group MI values to anatomical categories
        if Wh_theta
            eval(['MI_theta.' id1 '(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
        end
        if Wh_noth
            eval(['MI_noth.' id1 '(end+1:end+dsc) = fmh(DF_noth_eu,dsc);']);
        end
    end
    if ~isempty(id2)
        if Wh_theta
            eval(['MI_theta.' id2 '(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
        end
        if Wh_noth
            eval(['MI_noth.' id2 '(end+1:end+dsc) = fmh(DF_noth_eu,dsc);']);
        end
    end
    if isempty([id1 id2])
        if Wh_theta
            eval(['MI_theta.na(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
        end
        if Wh_noth
            eval(['MI_noth.na(end+1:end+dsc) = fmh(DF_noth_eu,dsc);']);
        end
    end
    
end

% Boxplot
H1 = mibox(MI_theta);
% H2 = mibox(MI_noth);

% Save
cd(resdir)
saveas(H1,'TEbox')
% saveas(H2,'MIbox_noth')
save('TEstat','MI_theta','MI_noth')
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
m1 = m1(~isnan(m1));
m2 = m2(~isnan(m2));
m3 = m3(~isnan(m3));
m4 = m4(~isnan(m4));
m5 = m5(~isnan(m5));

H = figure;
boxplot([m1 m2 m3 m4 m5],[zeros(size(m1)) ones(size(m2)) 2*ones(size(m3))...
    3*ones(size(m4)) 4*ones(size(m5))],'labels',[{'PV+'} {'PV-'} {'HCN+'}...
    {'HCN-'} {'na'}],'positions',[1 1.8 3.3 4.1 5.6]);
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

% X = X(X~=0);
X_ds = X(1:dsc);