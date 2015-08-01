function testat_burst
%TESTAT_BURST   Calculates statistics for transfer entropy data.
%   TESTAT_BURST runs on the output of TRTRENRUN. Edit code to specify
%   exact directories!
%
%   It runs on transfer entropy results, and addresses the question whether
%   transfer entropy is different among anatomically defined categories of
%   identified cells. It plots and saves a boxplot showing the transfer
%   entropy value distributions and the results of the statistical
%   comparisons by Mann-Whitney U-test.
%
%   TESTAT_BURST restricts the comparison to theta-bursting cells with MI
%   values significantly higher than control.
%
%   See also TESTAT and TESTAT_FPATTERN3.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_cont3\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_cont3\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_cont3\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_cont3\control\'];
clusdir_theta = [DATAPATH 'MI\Cluster\Theta\'];
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
DF_theta = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
    'HCNp',[],'HCNp_control',[],'HCNn',[],'HCNn_control',[],...
    'na',[],'na_control',[]);
NTE_thetaeu = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
    'HCNp',[],'HCNp_control',[],'HCNn',[],'HCNn_control',[],...
    'na',[],'na_control',[]);
NTE_thetaue = struct('PVp',[],'PVp_control',[],'PVn',[],'PVn_control',[],...
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
    
    if ~isequal(ftheta(1:6),fnoth(1:6))     % get identification data
        error('Technical error 40')
    else
        fname = ftheta(1:6);
    end
    inx = find(strcmp({tbl{:,1}},fname));
    id1 = tbl{inx,2};
    id2 = tbl{inx,3};
    
    cd(clusdir_theta)     % load burst data
    load([ftheta(1:end-15) 'CLUSTER'])
    if ClusterNumber == 0       % skip non-bursting cells
        continue
    end
    
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
    if ~isempty(id1)        % group MI values to anatomical categories
        if Wh_theta
            eval(['DF_theta.' id1 '(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
            eval(['NTE_thetaeu.' id1 '(end+1:end+dsc) = fmh(NTE_theta_eu,dsc);']);
            eval(['NTE_thetaue.' id1 '(end+1:end+dsc) = fmh(NTE_theta_ue,dsc);']);
        end
    end
    if ~isempty(id2)
        if Wh_theta
            eval(['DF_theta.' id2 '(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
            eval(['NTE_thetaeu.' id2 '(end+1:end+dsc) = fmh(NTE_theta_eu,dsc);']);
            eval(['NTE_thetaue.' id2 '(end+1:end+dsc) = fmh(NTE_theta_ue,dsc);']);
        end
    end
    if isempty([id1 id2])
        if Wh_theta
            eval(['DF_theta.na(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
            eval(['NTE_thetaeu.na(end+1:end+dsc) = fmh(NTE_theta_eu,dsc);']);
            eval(['NTE_thetaue.na(end+1:end+dsc) = fmh(NTE_theta_ue,dsc);']);
        end
    end
end

% Boxplot
H = mibox(DF_theta,NTE_thetaeu,NTE_thetaue);

% Save
cd(resdir)
saveas(H,'TEbox_burst_01')
% save('TEstat_burst','MI_theta')
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
function H = mibox(B,C,D)

m1 = B.PVp;
m2 = B.PVn;
m3 = B.HCNp;
m4 = B.HCNn;
m5 = B.na;
p1 = C.PVp;
p2 = C.PVn;
p3 = C.HCNp;
p4 = C.HCNn;
p5 = C.na;
q1 = D.PVp;
q2 = D.PVn;
q3 = D.HCNp;
q4 = D.HCNn;
q5 = D.na;
m1 = m1(~isnan(m1)&p1<=1&q1<=1);
m2 = m2(~isnan(m2)&p2<=1&q2<=1);
m3 = m3(~isnan(m3)&p3<=1&q3<=1);
m4 = m4(~isnan(m4)&p4<=1&q4<=1);
m5 = m5(~isnan(m5)&p5<=1&q5<=1);

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

X_ds = X(1:dsc);