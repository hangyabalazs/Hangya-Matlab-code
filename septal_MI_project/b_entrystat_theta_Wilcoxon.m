function [ThetaNoburst,ThetaBurst,TN,TB] = b_entrystat_theta_Wilcoxon
%ENTRYSTAT_THETA_WILCOXON   Calculates statistics for entropy data (theta segments).
%   ENTRYSTAT_THETA_WILCOXON runs on the output of ENTRYRUN_THETA or
%   ENTRYRUN3. Edit code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves Wilcoxon 
%   ranksum-test (Mann-Whitney U-test) results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_NOTH_WILCOXON.

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
    'W_hypothesis_eu',{},'W_significance_eu',{},'F_hypothesis_ue',{},...
    'F_significance_ue',{},'W_hypothesis_ue',{},'W_significance_ue',{});
ThetaBurst = struct('name',{},'F_hypothesis_eu',{},'F_significance_eu',{},...
    'W_hypothesis_eu',{},'W_significance_eu',{},'F_hypothesis_ue',{},...
    'F_significance_ue',{},'W_hypothesis_ue',{},'W_significance_ue',{});
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
    [Wp_eu,Wh_eu] = b_ranksum2(PDC_eu_real,PDC_eu_ctrl);
    [Wp_ue,Wh_ue] = b_ranksum2(PDC_ue_real,PDC_ue_ctrl);
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-8);
        ThetaNoburst(end).F_hypothesis_eu = Fh_eu;
        ThetaNoburst(end).F_significance_eu = Fp_eu;
        ThetaNoburst(end).W_hypothesis_eu = Wh_eu;
        ThetaNoburst(end).W_significance_eu = Wp_eu;
        ThetaNoburst(end).F_hypothesis_ue = Fh_ue;
        ThetaNoburst(end).F_significance_ue = Fp_ue;
        ThetaNoburst(end).W_hypothesis_ue = Wh_ue;
        ThetaNoburst(end).W_significance_ue = Wp_ue;
        ThetaNoburst(end).mean_eu_real = mPeur;
        ThetaNoburst(end).mean_eu_ctrl = mPeuc;
        ThetaNoburst(end).mean_ue_real = mPuer;
        ThetaNoburst(end).mean_ue_ctrl = mPuec;
        TN(end+1,1) = Fh_eu;
        TN(end,2) = Fp_eu;
        TN(end,3) = Wh_eu;
        TN(end,4) = Wp_eu;
        TN(end,5) = Fh_ue;
        TN(end,6) = Fp_ue;
        TN(end,7) = Wh_ue;
        TN(end,8) = Wp_ue;
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
        ThetaBurst(end).W_hypothesis_eu = Wh_eu;
        ThetaBurst(end).W_significance_eu = Wp_eu;
        ThetaBurst(end).F_hypothesis_ue = Fh_ue;
        ThetaBurst(end).F_significance_ue = Fp_ue;
        ThetaBurst(end).W_hypothesis_ue = Wh_ue;
        ThetaBurst(end).W_significance_ue = Wp_ue;
        ThetaBurst(end).mean_eu_real = mPeur;
        ThetaBurst(end).mean_eu_ctrl = mPeuc;
        ThetaBurst(end).mean_ue_real = mPuer;
        ThetaBurst(end).mean_ue_ctrl = mPuec;
        TB(end+1,1) = Fh_eu;
        TB(end,2) = Fp_eu;
        TB(end,3) = Wh_eu;
        TB(end,4) = Wp_eu;
        TB(end,5) = Fh_ue;
        TB(end,6) = Fp_ue;
        TB(end,7) = Wh_ue;
        TB(end,8) = Wp_ue;
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
fn = [DATAPATH 'Entry\Stat\theta_Wilcoxon_1sided.xls'];
str = [{'F_hypothesis_eu'} {'F_significance_eu'} {'W_hypothesis_eu'} {'W_significance_eu'}...
    {'F_hypothesis_ue'} {'F_significance_ue'} {'W_hypothesis_ue'} {'W_significance_ue'}...
    {'mean_PDC_eu_real'} {'mean_PDC_eu_ctrl'} {'eu_r>c'}...
    {'mean_PDC_ue_real'} {'mean_PDC_ue_ctrl'} {'ue_r>c'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');