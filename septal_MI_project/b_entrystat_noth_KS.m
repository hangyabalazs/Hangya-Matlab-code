function [NothNoburst,NothBurst,NN,NB] = b_entrystat_noth_KS
%ENTRYSTAT_NOTH_KS   Calculates statistics for entropy data (non-theta segments).
%   ENTRYSTAT_NOTH_KS runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact direcctories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves Kolmogorov-
%   Smirnov test results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_THETA_KS.

%Directories
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
NothNoburst = struct('name',{},'KS_hypothesis_eu',{},'KS_significance_eu',{},...
    'KS_hypothesis_ue',{},'KS_significance_ue',{});
NothBurst = struct('name',{},'KS_hypothesis_eu',{},'KS_significance_eu',{},...
    'KS_hypothesis_ue',{},'KS_significance_ue',{});
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
    [KSh_eu,KSp_eu] = kstest2(PDC_eu_real,PDC_eu_ctrl,0.05,'smaller');
    [KSh_ue,KSp_ue] = kstest2(PDC_ue_real,PDC_ue_ctrl,0.05,'smaller');
    if ClusterNumber == 0
        NothNoburst(end+1).name = files{o}(1:end-8);
        NothNoburst(end).KS_hypothesis_eu = KSh_eu;
        NothNoburst(end).KS_significance_eu = KSp_eu;
        NothNoburst(end).KS_hypothesis_ue = KSh_ue;
        NothNoburst(end).KS_significance_ue = KSp_ue;
        NothNoburst(end).mean_eu_real = mPeur;
        NothNoburst(end).mean_eu_ctrl = mPeuc;
        NothNoburst(end).mean_ue_real = mPuer;
        NothNoburst(end).mean_ue_ctrl = mPuec;
        NN(end+1,1) = KSh_eu;
        NN(end,2) = KSp_eu;
        NN(end,3) = KSh_ue;
        NN(end,4) = KSp_ue;
        NN(end,5) = mPeur;
        NN(end,6) = mPeuc;
        NN(end,7) = mPeur > mPeuc;
        NN(end,8) = mPuer;
        NN(end,9) = mPuec;
        NN(end,10) = mPuer > mPuec;
        NNname{end+1} = files{o}(1:end-8);
    else
        NothBurst(end+1).name = files{o}(1:end-8);
        NothBurst(end).KS_hypothesis_eu = KSh_eu;
        NothBurst(end).KS_significance_eu = KSp_eu;
        NothBurst(end).KS_hypothesis_ue = KSh_ue;
        NothBurst(end).KS_significance_ue = KSp_ue;
        NothBurst(end).mean_eu_real = mPeur;
        NothBurst(end).mean_eu_ctrl = mPeuc;
        NothBurst(end).mean_ue_real = mPuer;
        NothBurst(end).mean_ue_ctrl = mPuec;
        NB(end+1,1) = KSh_eu;
        NB(end,2) = KSp_eu;
        NB(end,3) = KSh_ue;
        NB(end,4) = KSp_ue;
        NB(end,5) = mPeur;
        NB(end,6) = mPeuc;
        NB(end,7) = mPeur > mPeuc;
        NB(end,8) = mPuer;
        NB(end,9) = mPuec;
        NB(end,10) = mPuer > mPuec;
        NBname{end+1} = files{o}(1:end-8);
    end
end

% Save
fn = [DATAPATH 'Entry\Stat\noth_KS_smaller.xls'];
str = [{'KS_hypothesis_eu'} {'KS_significance_eu'}...
    {'KS_hypothesis_ue'} {'KS_significance_ue'}...
    {'mean_PDC_eu_real'} {'mean_PDC_eu_ctrl'} {'eu_r>c'}...
    {'mean_PDC_ue_real'} {'mean_PDC_ue_ctrl'} {'ue_r>c'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');