function [Stats_eu, Stats_ue, STeu, STue] = emistat
%EMISTAT   Calculates statistics for mutual information data.
%   EMISTAT runs on the output of ENTRYRUN3C. Edit code to specify exact 
%   directories!
%
%   It runs on entropy results, and addresses the question whether mutual
%   information is different among anatomically defined categories of
%   identified cells. It plots and saves a boxplot showing the mutual value
%   distributions and the results of the statistical comparisons by
%   Mann-Whitney U-test.
%
%   See also ENTRYSTAT6.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta\control\'];
inpdir_theta_zshift = [DATAPATH 'MI\Theta_zshift\real\'];
inpdir_theta_control_zshift = [DATAPATH 'MI\Theta_zshift\control\'];
resdir = [DATAPATH 'MI\Stat\zshift\'];
mm = pwd;

% Files
[files_theta, files_short_theta] = genfiles(inpdir_theta);
files_short = files_short_theta;
sf = length(files_short);

% Load table
ff = [DATAPATH '5UEntropy\cells'];   % load identification data
[tbl0 tbl] = xlsread(ff);

% Statistics
MI_theta = [];
MI_theta_zshift = [];
for o = 1:sf
    cd(inpdir_theta)    % load MI values
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_real_theta = aIxy;
    
    cd(inpdir_theta_zshift)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_real_theta_zshift = aIxy;
    
    cd(inpdir_theta_control)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_control_theta = aIxy;
    
    cd(inpdir_theta_control_zshift)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_control_theta_zshift = aIxy;
    
    fname = ftheta(1:6);
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
    if median(MI_real_theta_zshift) > median(MI_control_theta_zshift)
        [Wp_theta_zshift,Wh_theta_zshift] = b_signrank2(MI_real_theta_zshift,MI_control_theta_zshift,'alpha',0.01);
    else
        Wh_theta_zshift = 0;
    end
    if (isequal(id1,'PVp') | isequal(id2,'HCNp')) && Wh_theta
        eval(['MI_theta(end+1:end+dsc) = fmh(MI_real_theta,dsc);']);
    end
    if (isequal(id1,'PVp') | isequal(id2,'HCNp')) && Wh_theta_zshift
        eval(['MI_theta_zshift(end+1:end+dsc) = fmh(MI_real_theta_zshift,dsc);']);
    end
end

% Boxplot
H = mibox(MI_theta,MI_theta_zshift,'original PV or HCN','shifted PV or HCN');

% Save
cd(resdir)
saveas(H,'MIbox_zshift_PVHCN')
save('MIstat_zshift_PVHCN','MI_theta','MI_theta_zshift')
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