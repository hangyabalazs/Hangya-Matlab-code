function [ThetaNoburst,ThetaBurst,TN,TB] = b_entrystat_theta_Ft
%ENTRYSTAT_NOTH_FT   Calculates statistics for entropy data (theta segments).
%   ENTRYSTAT_NOTH_FT runs on the output of ENTRYRUN_THETA or ENTRYRUN3.
%   Edit code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves F- and t-test
%   results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_NOTH_T.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry\'];
clusdir = [DATAPATH 'Burst\Cluster\Theta\'];
mm = pwd;
cd(inpdir)
cd(['Theta\dtf\halfsec\real'])
files0 = dir(pwd);
files = {};
files_clus = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if ext == '.mat'
        files{end+1} = files0(i).name;
        files_clus{end+1} = [files0(i).name(1:end-7) 'CLUSTER'];
    end
end
sf = length(files);

% Statistics
ThetaNoburst = struct('name',{},'F_hypothesis_eu',{},'F_significance_eu',{},...
    't_hypothesis_eu',{},'t_significance_eu',{},'F_hypothesis_ue',{},...
    'F_significance_ue',{},'t_hypothesis_ue',{},'t_significance_ue',{});
ThetaBurst = struct('name',{},'F_hypothesis_eu',{},'F_significance_eu',{},...
    't_hypothesis_eu',{},'t_significance_eu',{},'F_hypothesis_ue',{},...
    'F_significance_ue',{},'t_hypothesis_ue',{},'t_significance_ue',{});
TN = [];
TB = [];
TNname = {};
TBname = {};
for o = 1:sf
    cd(inpdir)
    cd(['Theta\dtf\halfsec\real'])
    load(files{o})
    PDC_eu_real = mean(PDC_eu(3:6,:));
    PDC_ue_real = mean(PDC_ue(3:6,:));
    mPeur = mean(PDC_eu_real);
    mPuer = mean(PDC_ue_real);
    mm2 = pwd;
    cd(clusdir)
    load(files_clus{o})
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    PDC_eu_ctrl = mean(PDC_eu(3:6,:));
    PDC_ue_ctrl = mean(PDC_ue(3:6,:));
    mPeuc = mean(PDC_eu_ctrl);
    mPuec = mean(PDC_ue_ctrl);
    [stat_eu,Fh_eu,Fp_eu] = b_Ftest(PDC_eu_real,PDC_eu_ctrl);
    [stat_ue,Fh_ue,Fp_ue] = b_Ftest(PDC_ue_real,PDC_ue_ctrl);
    [th_eu,tp_eu,ci_eu] = ttest2(PDC_eu_real,PDC_eu_ctrl,[],'right','unequal');
    [th_ue,tp_ue,ci_ue] = ttest2(PDC_ue_real,PDC_ue_ctrl,[],'right','unequal');
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-8);
        ThetaNoburst(end).F_hypothesis_eu = Fh_eu;
        ThetaNoburst(end).F_significance_eu = Fp_eu;
        ThetaNoburst(end).t_hypothesis_eu = th_eu;
        ThetaNoburst(end).t_significance_eu = tp_eu;
        ThetaNoburst(end).F_hypothesis_ue = Fh_ue;
        ThetaNoburst(end).F_significance_ue = Fp_ue;
        ThetaNoburst(end).t_hypothesis_ue = th_ue;
        ThetaNoburst(end).t_significance_ue = tp_ue;
        ThetaNoburst(end).mean_eu_real = mPeur;
        ThetaNoburst(end).mean_eu_ctrl = mPeuc;
        ThetaNoburst(end).mean_ue_real = mPuer;
        ThetaNoburst(end).mean_ue_ctrl = mPuec;
        TN(end+1,1) = Fh_eu;
        TN(end,2) = Fp_eu;
        TN(end,3) = th_eu;
        TN(end,4) = tp_eu;
        TN(end,5) = Fh_ue;
        TN(end,6) = Fp_ue;
        TN(end,7) = th_ue;
        TN(end,8) = tp_ue;
        TN(end,9) = mPeur;
        TN(end,10) = mPeuc;
        TN(end,11) = mPeur > mPeuc;
        TN(end,12) = mPuer;
        TN(end,13) = mPuec;
        TN(end,14) = mPuer > mPuec;
        TNname{end+1} = files{o}(1:end-8);
    else
        ThetaBurst(end+1).name = files{o}(1:end-8);
        ThetaBurst(end).F_hypothesis_eu = Fh_eu;
        ThetaBurst(end).F_significance_eu = Fp_eu;
        ThetaBurst(end).t_hypothesis_eu = th_eu;
        ThetaBurst(end).t_significance_eu = tp_eu;
        ThetaBurst(end).F_hypothesis_ue = Fh_ue;
        ThetaBurst(end).F_significance_ue = Fp_ue;
        ThetaBurst(end).t_hypothesis_ue = th_ue;
        ThetaBurst(end).t_significance_ue = tp_ue;
        ThetaBurst(end).mean_eu_real = mPeur;
        ThetaBurst(end).mean_eu_ctrl = mPeuc;
        ThetaBurst(end).mean_ue_real = mPuer;
        ThetaBurst(end).mean_ue_ctrl = mPuec;
        TB(end+1,1) = Fh_eu;
        TB(end,2) = Fp_eu;
        TB(end,3) = th_eu;
        TB(end,4) = tp_eu;
        TB(end,5) = Fh_ue;
        TB(end,6) = Fp_ue;
        TB(end,7) = th_ue;
        TB(end,8) = tp_ue;
        TB(end,9) = mPeur;
        TB(end,10) = mPeuc;
        TB(end,11) = mPeur > mPeuc;
        TB(end,12) = mPuer;
        TB(end,13) = mPuec;
        TB(end,14) = mPuer > mPuec;
        TBname{end+1} = files{o}(1:end-8);
    end
end

% Save
fn = [DATAPATH 'Entry\Stat\theta_tailed.xls'];
str = [{'F_hypothesis_eu'} {'F_significance_eu'} {'t_hypothesis_eu'} {'t_significance_eu'}...
    {'F_hypothesis_ue'} {'F_significance_ue'} {'t_hypothesis_ue'} {'t_significance_ue'}...
    {'mean_PDC_eu_real'} {'mean_PDC_eu_ctrl'} {'eu_r>c'}...
    {'mean_PDC_ue_real'} {'mean_PDC_ue_ctrl'} {'ue_r>c'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');