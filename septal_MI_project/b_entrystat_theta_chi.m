function [ThetaNoburst,ThetaBurst,TN,TB] = b_entrystat_theta_chi
%ENTRYSTAT_THETA_CHI   Calculates statistics for entropy data (theta segments).
%   ENTRYSTAT_THETA_CHI runs on the output of ENTRYRUN_THETA or ENTRYRUN3.
%   Edit code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves chi-square-test 
%   results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_THETA_FT.

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
ThetaNoburst = struct('name',{},'chi_hypothesis_eu',{},'chi_significance_eu',{},...
    'KS_hypothesis_eu',{},'KS_significance_eu',{},'chi_hypothesis_ue',{},...
    'chi_significance_ue',{},'KS_hypothesis_ue',{},'KS_significance_ue',{});
ThetaBurst = struct('name',{},'chi_hypothesis_eu',{},'chi_significance_eu',{},...
    'KS_hypothesis_eu',{},'KS_significance_eu',{},'chi_hypothesis_ue',{},...
    'chi_significance_ue',{},'KS_hypothesis_ue',{},'KS_significance_ue',{});
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
    [ch_eu,cp_eu] = b_chi2test(PDC_eu_real,PDC_eu_ctrl);
    [ch_ue,cp_ue] = b_chi2test(PDC_ue_real,PDC_ue_ctrl);
    [KSh_eu,KSp_eu] = kstest2(PDC_eu_real,PDC_eu_ctrl,0.05,'larger');
    [KSh_ue,KSp_ue] = kstest2(PDC_ue_real,PDC_ue_ctrl,0.05,'larger');
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-8);
        ThetaNoburst(end).chi_hypothesis_eu = ch_eu;
        ThetaNoburst(end).chi_significance_eu = cp_eu;
        ThetaNoburst(end).KS_hypothesis_eu = KSh_eu;
        ThetaNoburst(end).KS_significance_eu = KSp_eu;
        ThetaNoburst(end).chi_hypothesis_ue = ch_ue;
        ThetaNoburst(end).chi_significance_ue = cp_ue;
        ThetaNoburst(end).KS_hypothesis_ue = KSh_ue;
        ThetaNoburst(end).KS_significance_ue = KSp_ue;
        ThetaNoburst(end).mean_eu_real = mPeur;
        ThetaNoburst(end).mean_eu_ctrl = mPeuc;
        ThetaNoburst(end).mean_ue_real = mPuer;
        ThetaNoburst(end).mean_ue_ctrl = mPuec;
        TN(end+1,1) = ch_eu;
        TN(end,2) = cp_eu;
        TN(end,3) = KSh_eu;
        TN(end,4) = KSp_eu;
        TN(end,5) = ch_ue;
        TN(end,6) = cp_ue;
        TN(end,7) = KSh_ue;
        TN(end,8) = KSp_ue;
        TN(end,9) = mPeur;
        TN(end,10) = mPeuc;
        TN(end,11) = mPeur > mPeuc;
        TN(end,12) = mPuer;
        TN(end,13) = mPuec;
        TN(end,14) = mPuer > mPuec;
        TNname{end+1} = files{o}(1:end-8);
    else
        ThetaBurst(end+1).name = files{o}(1:end-8);
        ThetaBurst(end).chi_hypothesis_eu = ch_eu;
        ThetaBurst(end).chi_significance_eu = cp_eu;
        ThetaBurst(end).KS_hypothesis_eu = KSh_eu;
        ThetaBurst(end).KS_significance_eu = KSp_eu;
        ThetaBurst(end).chi_hypothesis_ue = ch_ue;
        ThetaBurst(end).chi_significance_ue = cp_ue;
        ThetaBurst(end).KS_hypothesis_ue = KSh_ue;
        ThetaBurst(end).KS_significance_ue = KSp_ue;
        ThetaBurst(end).mean_eu_real = mPeur;
        ThetaBurst(end).mean_eu_ctrl = mPeuc;
        ThetaBurst(end).mean_ue_real = mPuer;
        ThetaBurst(end).mean_ue_ctrl = mPuec;
        TB(end+1,1) = ch_eu;
        TB(end,2) = cp_eu;
        TB(end,3) = KSh_eu;
        TB(end,4) = KSp_eu;
        TB(end,5) = ch_ue;
        TB(end,6) = cp_ue;
        TB(end,7) = KSh_ue;
        TB(end,8) = KSp_ue;
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
fn = [DATAPATH 'Entry\Stat\theta_chiKS.xls'];
str = [{'chi_hypothesis_eu'} {'chi_significance_eu'} {'KS_hypothesis_eu'} {'KS_significance_eu'}...
    {'chi_hypothesis_ue'} {'chi_significance_ue'} {'KS_hypothesis_ue'} {'KS_significance_ue'}...
    {'mean_PDC_eu_real'} {'mean_PDC_eu_ctrl'} {'eu_r>c'}...
    {'mean_PDC_ue_real'} {'mean_PDC_ue_ctrl'} {'ue_r>c'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');