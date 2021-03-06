function [NothNoburst,NothBurst,NN,NB] = b_entrystat_noth_Wilcoxonmod_1sided
%ENTRYSTAT_NOTH_WILCOXONMOD_1SIDED   Calculates statistics for entropy data (non-theta segments).
%   ENTRYSTAT_NOTH_WILCOXONMOD_1SIDED runs on the output of ENTRYRUN_THETA
%   or ENTRYRUN3. Edit code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves one-sided
%   Wilcoxon ranksum-test (Mann-Whitney U-test) results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_THETA_WILCOXONMOD_1SIDED.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry\'];
clusdir = [DATAPATH 'Burst\Cluster\Noth\'];
mm = pwd;
cd(inpdir)
cd(['Noth\dtf\halfsec\real'])
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
NothNoburst = struct('name',{},'F_hypothesis_eu',{},'F_significance_eu',{},...
    'W_hypothesis_eu',{},'W_significance_eu',{},'F_hypothesis_ue',{},...
    'F_significance_ue',{},'W_hypothesis_ue',{},'W_significance_ue',{});
NothBurst = struct('name',{},'F_hypothesis_eu',{},'F_significance_eu',{},...
    'W_hypothesis_eu',{},'W_significance_eu',{},'F_hypothesis_ue',{},...
    'F_significance_ue',{},'W_hypothesis_ue',{},'W_significance_ue',{});
NN = [];
NB = [];
NNname = {};
NBname = {};
for o = 1:sf
    cd(inpdir)
    cd(['Noth\dtf\halfsec\real'])
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
    if mPeuc > mPeur
        Wh_eu = 0;
        Wp_eu = -1;
    end
    if mPuec > mPuer
        Wh_ue = 0;
        Wp_ue = -1;
    end
    if ClusterNumber == 0
        NothNoburst(end+1).name = files{o}(1:end-8);
        NothNoburst(end).F_hypothesis_eu = Fh_eu;
        NothNoburst(end).F_significance_eu = Fp_eu;
        NothNoburst(end).W_hypothesis_eu = Wh_eu;
        NothNoburst(end).W_significance_eu = Wp_eu;
        NothNoburst(end).F_hypothesis_ue = Fh_ue;
        NothNoburst(end).F_significance_ue = Fp_ue;
        NothNoburst(end).W_hypothesis_ue = Wh_ue;
        NothNoburst(end).W_significance_ue = Wp_ue;
        NothNoburst(end).mean_eu_real = mPeur;
        NothNoburst(end).mean_eu_ctrl = mPeuc;
        NothNoburst(end).mean_ue_real = mPuer;
        NothNoburst(end).mean_ue_ctrl = mPuec;
        NN(end+1,1) = Fh_eu;
        NN(end,2) = Fp_eu;
        NN(end,3) = Wh_eu;
        NN(end,4) = Wp_eu;
        NN(end,5) = Fh_ue;
        NN(end,6) = Fp_ue;
        NN(end,7) = Wh_ue;
        NN(end,8) = Wp_ue;
        NN(end,9) = mPeur;
        NN(end,10) = mPeuc;
        NN(end,11) = mPeur > mPeuc;
        NN(end,12) = mPuer;
        NN(end,13) = mPuec;
        NN(end,14) = mPuer > mPuec;
        NNname{end+1} = files{o}(1:end-8);
    else
        NothBurst(end+1).name = files{o}(1:end-8);
        NothBurst(end).F_hypothesis_eu = Fh_eu;
        NothBurst(end).F_significance_eu = Fp_eu;
        NothBurst(end).W_hypothesis_eu = Wh_eu;
        NothBurst(end).W_significance_eu = Wp_eu;
        NothBurst(end).F_hypothesis_ue = Fh_ue;
        NothBurst(end).F_significance_ue = Fp_ue;
        NothBurst(end).W_hypothesis_ue = Wh_ue;
        NothBurst(end).W_significance_ue = Wp_ue;
        NothBurst(end).mean_eu_real = mPeur;
        NothBurst(end).mean_eu_ctrl = mPeuc;
        NothBurst(end).mean_ue_real = mPuer;
        NothBurst(end).mean_ue_ctrl = mPuec;
        NB(end+1,1) = Fh_eu;
        NB(end,2) = Fp_eu;
        NB(end,3) = Wh_eu;
        NB(end,4) = Wp_eu;
        NB(end,5) = Fh_ue;
        NB(end,6) = Fp_ue;
        NB(end,7) = Wh_ue;
        NB(end,8) = Wp_ue;
        NB(end,9) = mPeur;
        NB(end,10) = mPeuc;
        NB(end,11) = mPeur > mPeuc;
        NB(end,12) = mPuer;
        NB(end,13) = mPuec;
        NB(end,14) = mPuer > mPuec;
        NBname{end+1} = files{o}(1:end-8);
    end
end

% Save
fn = [DATAPATH 'Entry\Stat\noth_Wilcoxonmod_1sided.xls'];
str = [{'F_hypothesis_eu'} {'F_significance_eu'} {'W_hypothesis_eu'} {'W_significance_eu'}...
    {'F_hypothesis_ue'} {'F_significance_ue'} {'W_hypothesis_ue'} {'W_significance_ue'}...
    {'mean_PDC_eu_real'} {'mean_PDC_eu_ctrl'} {'eu_r>c'}...
    {'mean_PDC_ue_real'} {'mean_PDC_ue_ctrl'} {'ue_r>c'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');