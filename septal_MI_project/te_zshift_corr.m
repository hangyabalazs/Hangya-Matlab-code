function te_zshift_corr
%TE_ZSHIFT_CORR   Calculates relation of transfer entropy and zshift.
%   TE_ZSHIFT_CORR runs on the output of TRTRENRUN and EZSHIFTRUN. Edit
%   code to specify exact directories!
%
%   It runs on transfer entropy results, and addresses the question whether
%   transfer entropy is correlated with zshift in anatomically defined
%   categories of identified cells.
%
%   See also TESTAT, TRTRENRUN and EZSHIFTRUN.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'MI\Theta_cont3\real\'];
inpdir_theta_control = [DATAPATH 'MI\Theta_cont3\control\'];
inpdir_te = [DATAPATH 'TE\timeresolved_tauNTEmax\'];
resdir = [DATAPATH 'TE\Stat_cont3_tauNTEmax\'];
mm = pwd;

% Files
[files_theta, files_short] = genfiles(inpdir_theta);
sf = length(files_short);

% Load table
ff = [DATAPATH '5UEntropy\cells'];   % load identification data
[tbl0 tbl] = xlsread(ff);

% Read Z-shift excel file
dr = [DATAPATH 'Ezshift\'];
fn = [dr 'zshift_longest_theta.xls'];
cd(dr)
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'maxlen');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Statistics
H1 = figure;
hold on
H2 = figure;
hold on
PVp_zsh = [];
PVn_zsh = [];
PVna_zsh = [];
HCNp_zsh = [];
HCNn_zsh = [];
PVp_MI = [];
PVn_MI = [];
PVna_MI = [];
HCNp_MI = [];
HCNn_MI = [];
for o = 1:sf
    cd(inpdir_theta)    % load MI values
    ftheta = files_theta{o};
    load(ftheta)
    MI_real_theta = aIxy;
        
    cd(inpdir_theta_control)
    inx = find(strcmp(files_short{o},files_short));
    ftheta = files_theta{inx};
    load(ftheta)
    MI_control_theta = aIxy;
    
    fname = ftheta(1:6);    % find anatomical identification of the cells
    inx = find(strcmp(fname,tbl(:,1)));
    pv = tbl{inx,2};
    hcn = tbl{inx,3};
    
    inx = find(strcmp(fname,atx(:,1)));     % find Z-shift
    zsh = atx{inx,4};
    
    fs = findstr(ftheta,'_');       % load transfer entropy
    ff = [inpdir_te ftheta(1:fs(2)-1) '_TE.mat'];
    load(ff)
    mMI = median(DF_theta_eu(1:4));
    
    dsc = 4;
    fmh = @select;
    if median(MI_real_theta) > median(MI_control_theta)     % check if MI > control
        [Wp_theta,Wh_theta] = b_signrank2(MI_real_theta,MI_control_theta,'alpha',0.01);
    else
        Wh_theta = 0;
    end
    
    figure(H1)      % plot
    if Wh_theta && ~isnan(zsh)
        if isequal(pv,'PVp')
            h1 = plot(zsh,mMI,'r.','MarkerSize',50);
            PVp_zsh(end+1) = zsh;
            PVp_MI(end+1) = mMI;
        elseif isequal(pv,'PVn')
            h2 = plot(zsh,mMI,'.','MarkerSize',50,'MarkerFaceColor',...
                [1 0.69 0.39],'MarkerEdgeColor',[1 0.69 0.39]);
            PVn_zsh(end+1) = zsh;
            PVn_MI(end+1) = mMI;
        elseif isempty(pv)
            h3 = plot(zsh,mMI,'k.','MarkerSize',30);
            PVna_zsh(end+1) = zsh;
            PVna_MI(end+1) = mMI;
        end
    end
    figure(H2)
    if Wh_theta && ~isnan(zsh)
        if isequal(hcn,'HCNp')
            h1 = plot(zsh,mMI,'b.','MarkerSize',50);
            HCNp_zsh(end+1) = zsh;
            HCNp_MI(end+1) = mMI;
        elseif isequal(hcn,'HCNn')
            h2 = plot(zsh,mMI,'.','MarkerSize',50,'MarkerFaceColor',...
            [51 204 255]/256,'MarkerEdgeColor',[51 204 255]/256);
            HCNn_zsh(end+1) = zsh;
            HCNn_MI(end+1) = mMI;
        elseif isempty(hcn)
            h3 = plot(zsh,mMI,'k.','MarkerSize',30);
        end
    end
end

x = PVp_zsh';
y = PVp_MI';
X = [ones(size(x)) x];
[b,bint,r,rint,stats] = regress(y,X);
corrcoef(x,y);
R = sqrt(stats(1))         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3)

x = HCNp_zsh';
y = HCNp_MI';
X = [ones(size(x)) x];
[b,bint,r,rint,stats] = regress(y,X);
corrcoef(x,y);
R = sqrt(stats(1))         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3)

x = [PVp_zsh'; PVn_zsh'; PVna_zsh'];
y = [PVp_MI'; PVn_MI'; PVna_MI'];
X = [ones(size(x)) x];
[b,bint,r,rint,stats] = regress(y,X);
corrcoef(x,y);
R = sqrt(stats(1))         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3)

% Save
cd(resdir)
saveas(H1,'zshift_vs_MI_PV')
saveas(H2,'zshift_vs_MI_HCN')
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

X_ds = X(1:4);