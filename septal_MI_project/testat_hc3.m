function testat_hc3
%TESTAT_HC3   Calculates statistics for transfer entropy data.
%   TESTAT_HC3 runs on the output of TRTRENRUN. Edit code to specify exact 
%   directories!
%
%   It runs on entropy results, and addresses the question whether transfer
%   entropy is different among anatomically defined categories of
%   identified cells. It plots and saves a boxplot showing the transfer
%   entropy value distributions and the results of the statistical
%   comparisons by Mann-Whitney U-test.
%
%   See also TESTAT_FPATTERN3_HC and TRTRENRUN.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_hc\real\'];
inpdir_noth = [DATAPATH 'MI\Noth_hc\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_hc\control\'];
inpdir_noth_control = [DATAPATH 'MI\Noth_hc\control\'];
inpdir_te = [DATAPATH 'TE\timeresolved_tauNTEmax_hc\'];
resdir = [DATAPATH 'TE\Stat_cont3_tauNTEmax_hc\'];
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
    
     fs = findstr(ftheta,'_');       % load transfer entropy
%     ff = [inpdir_te ftheta(1:fs(2)-1) '_TE.mat'];
    ff = [inpdir_te ftheta(1:6) '_TE.mat'];
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
            if Wp
                counter_noth = counter_noth + 1;
            end
        end
    end
    if Wh
        eval(['MI_theta(end+1:end+dsc) = fmh(DF_theta_eu,dsc);']);
        eval(['MI_noth(end+1:end+dsc) = fmh(DF_noth_eu,dsc);']);
    end
end

% Boxplot
H = mibox(MI_theta,MI_noth,'HC theta','HC noth');

% Save
cd(resdir)
saveas(H,'TEsegments3_hc')
save('TEstat3','MI_theta','MI_noth')
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

m1 = m1(~isnan(m1));
m2 = m2(~isnan(m2));
H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{l1} {l2}],...
    'positions',[0.5 1.1]);
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

X_ds = X(1:dsc);