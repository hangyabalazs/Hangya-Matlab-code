function [Theta,T] = b_2entry3chstat3_theta
%2ENTRY3CHSTAT3_THETA   Calculates statistics for entropy data (theta segments).
%   2ENTRY3CHSTAT3_THETA runs on the output of ENTRYRUN_3CH. Edit code to
%   specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   This program is for unit-unit analysis!
%
%   See also 2ENTRY3CHSTAT3_NOTH and ENTRYSTAT3_THETA.

% DTF
% theta
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry_uu_stand\'];
inpplus = ['Theta\dtf\onesec\real\'];
fn = [DATAPATH '2Entry_uu_stand\Stat\DTF_onesec\2theta.xls'];
fnm = [DATAPATH '2Entry_uu_stand\Stat\DTF_onesec\stat_2theta.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Theta_3ch\'];
mm = pwd;
cd(inpdir)
cd(inpplus)
files0 = dir(pwd);
files = {};
files_clus = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if isequal(ext,'.mat')
        files{end+1} = files0(i).name;
        files_clus{end+1} = [files0(i).name(1:end-7) 'CLUSTER'];
    end
end
sf = length(files);

% Statistics
Theta = struct([]);
T = {};
Tname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    PDC_eu_real = mean(PDC_eu(3:6,:));
    PDC_ue_real = mean(PDC_ue(3:6,:));
    mPeur = mean(PDC_eu_real);
    mPuer = mean(PDC_ue_real);
    mm2 = pwd;
    cd(clusdir)
    clusH = [files_clus{o}(1:3) 'H' files_clus{o}(5:end)];
    load(clusH)
    ClusterH = ClusterNumber;
    clusM = [files_clus{o}(1:3) 'M' files_clus{o}(5:end)];
    load(clusM)
    ClusterM = ClusterNumber;
    clear ClusterNumber
    cd(mm2)
    if mPeur > mPuer
        [KSh,KSp] = kstest2(PDC_eu_real,PDC_ue_real,0.05,'smaller');
        [th,tp,ci] = ttest(PDC_eu_real-PDC_ue_real,0,[],'right');
        [Wp,Wh] = b_signrank2(PDC_eu_real-PDC_ue_real);
    else
        [KSh,KSp] = kstest2(PDC_eu_real,PDC_ue_real,0.05,'larger');
        [th,tp,ci] = ttest(PDC_eu_real-PDC_ue_real,0,[],'left');
        [Wp,Wh] = b_signrank2(PDC_ue_real-PDC_eu_real);
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    Theta(end+1).name = files{o}(1:end-8);
    Theta(end).t_hypothesis = th;
    Theta(end).t_significance = tp;
    Theta(end).W_hypothesis = Wh;
    Theta(end).W_significance = Wp;
    Theta(end).KS_hypothesis = KSh;
    Theta(end).KS_significance = KSp;
    Theta(end).mean_eu_real = mPeur;
    Theta(end).mean_ue_real = mPuer;
    Theta(end).HC_unit_cluternumber = ClusterH;
    Theta(end).MS_unit_cluternumber = ClusterM;
    T{end+1,1} = th;
    T{end,2} = Wh;
    T{end,3} = KSh;
    T{end,4} = [];
    T{end,5} = mPeur > mPuer;
    T{end,5} = T{end,5} + 1 - 1;
    T{end,6} = [];
    T{end,7} = ClusterH;
    T{end,8} = ClusterM;
    Tname{end+1} = files{o}(1:end-8);
end

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'eu>?ue'} {' '} {'HC_cluster'} {'MS_cluster'}];
xlswrite(fn,str,'Theta','B1');
xlswrite(fn,Tname','Theta','A2');
xlswrite(fn,T,'Theta','B2');

save(fnm,'Theta','T')