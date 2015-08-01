function [ThetaNoburst,ThetaBurst,TN,TB] = b_2entrymmstat3_theta
%2ENTRYMMSTAT3_THETA   Calculates statistics for entropy data (theta segments).
%   2ENTRYMMSTAT3_THETA runs on the output of ENTRYRUN_THETA or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   It compares real ue-eu difference to control difference.
%
%   See also 2ENTRYSTAT3_NOTH and ENTRYSTAT4_THETA.

% DTF
% theta
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3b\'];
inpplus = ['Theta\dtf\onesec\real\'];
fn = [DATAPATH 'Entry3b_GG\Stat\DTF_onesec\2theta_ue.xls'];
fnm = [DATAPATH 'Entry3b_GG\Stat\DTF_onesec\stat_2theta_ue.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Theta\'];
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

%Statistics
ThetaNoburst = struct('name',{},'t_hypothesis_diff',{},'W_hypothesis_diff',{},...
    'KS_hypothesis_diff',{});
ThetaBurst = struct('name',{},'t_hypothesis_diff',{},'W_hypothesis_diff',{},...
    'KS_hypothesis_diff',{});
TN = {};
TB = {};
TNname = {};
TBname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    PDC_eu_real = mean(PDC_eu(3:6,:));
    PDC_ue_real = mean(PDC_ue(3:6,:));
    PDC_diff_real = PDC_ue_real - PDC_eu_real;
    mPdiffr = mean(PDC_diff_real);
    mm2 = pwd;
    cd(clusdir)
    load(files_clus{o})
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    PDC_eu_ctrl = mean(PDC_eu(3:6,:));
    PDC_ue_ctrl = mean(PDC_ue(3:6,:));
    PDC_diff_ctrl = PDC_ue_ctrl - PDC_eu_ctrl;
    mPdiffc = mean(PDC_diff_ctrl);
    alfa = 0.007;
    [KSh_diff,KSp_diff] = kstest2(PDC_diff_real,PDC_diff_ctrl,alfa,'smaller');
    [th_diff,tp_diff,ci_diff] = ttest(PDC_diff_real-PDC_diff_ctrl,0,alfa,'right');
    [Wp_diff,Wh_diff] = b_signrank2(PDC_diff_real,PDC_diff_ctrl,'alpha',alfa);
    KSh_diff = KSh_diff + 1 - 1;    % convert to numbers from logical values
    th_diff = th_diff + 1 - 1;
    Wh_diff = Wh_diff + 1 - 1;
    if mPdiffc > mPdiffr
        Wh_diff = 0;
        Wp_diff = -1;
    end
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-8);
        ThetaNoburst(end).t_hypothesis_diff = th_diff;
        ThetaNoburst(end).t_significance_diff = tp_diff;
        ThetaNoburst(end).W_hypothesis_diff = Wh_diff;
        ThetaNoburst(end).W_significance_diff = Wp_diff;
        ThetaNoburst(end).KS_hypothesis_diff = KSh_diff;
        ThetaNoburst(end).KS_significance_diff = KSp_diff;
        ThetaNoburst(end).mean_diff_real = mPdiffr;
        ThetaNoburst(end).mean_diff_ctrl = mPdiffc;
        TN{end+1,1} = th_diff;
        TN{end,2} = Wh_diff;
        TN{end,3} = KSh_diff;
        TN{end,4} = [];
        TN{end,5} = mPdiffr > mPdiffc;
        TN{end,5} = TN{end,5} + 1 - 1;
        TNname{end+1} = files{o}(1:end-8);
    else
        ThetaBurst(end+1).name = files{o}(1:end-8);
        ThetaBurst(end).t_hypothesis_diff = th_diff;
        ThetaBurst(end).t_significance_diff = tp_diff;
        ThetaBurst(end).W_hypothesis_diff = Wh_diff;
        ThetaBurst(end).W_significance_diff = Wp_diff;
        ThetaBurst(end).KS_hypothesis_diff = KSh_diff;
        ThetaBurst(end).KS_significance_diff = KSp_diff;
        ThetaBurst(end).mean_diff_real = mPdiffr;
        ThetaBurst(end).mean_diff_ctrl = mPdiffc;
        TB{end+1,1} = th_diff;
        TB{end,2} = Wh_diff;
        TB{end,3} = KSh_diff;
        TB{end,4} = [];
        TB{end,5} = mPdiffr > mPdiffc;
        TB{end,5} = TB{end,5} + 1 - 1;
        TBname{end+1} = files{o}(1:end-8);
    end
end
alfa

% Save
str = [{'t_hypothesis_diff'} {'W_hypothesis_diff'} {'KS_hypothesis_diff'}...
    {' '} {'diff_r>c'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');

save(fnm,'ThetaNoburst','ThetaBurst','TN','TB')
cd(mm)