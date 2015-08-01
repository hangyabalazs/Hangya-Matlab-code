function [Noth,N] = b_2entry3chstat_noth
%2ENTRY3CHSTAT_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRY3CHSTAT_NOTH runs on the output of ENTRYRUN_3CH. Edit code to
%   specify exact direcctories!
%
%   It runs on PDC results, and addresses the question whether real data is
%   larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   This program is for unit-unit analysis!
%
%   See also 2ENTRY3CHSTAT_THETA and ENTRYSTAT_NOTH.

% DTF
% noth
% eu r>c and ue r>c

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3c_unitunit\'];
inpplus = ['Noth\dtf\onesec\real'];
fn = [DATAPATH 'Entry3c_unitunit\Stat\DTF_onesec\noth_rorc.xls'];
fnm = [DATAPATH 'Entry3c_unitunit\Stat\DTF_onesec\stat_noth_rorc.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Noth_3ch\'];
mm = pwd;
cd(inpdir)
cd(inpplus)
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
Noth = struct([]);
N = {};
Nname = {};
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
    cd ..
    cd(['control'])
    load(files{o})
    PDC_eu_ctrl = mean(PDC_eu(3:6,:));
    PDC_ue_ctrl = mean(PDC_ue(3:6,:));
    mPeuc = mean(PDC_eu_ctrl);
    mPuec = mean(PDC_ue_ctrl);
    [KSh_eu,KSp_eu] = kstest2(PDC_eu_real,PDC_eu_ctrl,0.05,'smaller');
    [KSh_ue,KSp_ue] = kstest2(PDC_ue_real,PDC_ue_ctrl,0.05,'smaller');
    [th_eu,tp_eu,ci_eu] = ttest(PDC_eu_real-PDC_eu_ctrl,0,[],'right');
    [th_ue,tp_ue,ci_ue] = ttest(PDC_ue_real-PDC_ue_ctrl,0,[],'right');
    [Wp_eu,Wh_eu] = b_signrank2(PDC_eu_real-PDC_eu_ctrl);
    [Wp_ue,Wh_ue] = b_signrank2(PDC_ue_real-PDC_ue_ctrl);
    KSh_eu = KSh_eu + 1 - 1;    % convert to numbers from logical values
    th_eu = th_eu + 1 - 1;
    Wh_eu = Wh_eu + 1 - 1;
    KSh_ue = KSh_ue + 1 - 1;
    th_ue = th_ue + 1 - 1;
    Wh_ue = Wh_ue + 1 - 1;
    if mPeuc > mPeur
        Wh_eu = 0;
        Wp_eu = -1;
    end
    if mPuec > mPuer
        Wh_ue = 0;
        Wp_ue = -1;
    end
    Noth(end+1).name = files{o}(1:end-8);
    Noth(end).t_hypothesis_eu = th_eu;
    Noth(end).t_significance_eu = tp_eu;
    Noth(end).W_hypothesis_eu = Wh_eu;
    Noth(end).W_significance_eu = Wp_eu;
    Noth(end).KS_hypothesis_eu = KSh_eu;
    Noth(end).KS_significance_eu = KSp_eu;
    Noth(end).t_hypothesis_ue = th_ue;
    Noth(end).t_significance_ue = tp_ue;
    Noth(end).W_hypothesis_ue = Wh_ue;
    Noth(end).W_significance_ue = Wp_ue;
    Noth(end).KS_hypothesis_ue = KSh_ue;
    Noth(end).KS_significance_ue = KSp_ue;
    Noth(end).mean_eu_real = mPeur;
    Noth(end).mean_eu_ctrl = mPeuc;
    Noth(end).mean_ue_real = mPuer;
    Noth(end).mean_ue_ctrl = mPuec;
    N{end+1,1} = th_eu;
    N{end,2} = Wh_eu;
    N{end,3} = KSh_eu;
    N{end,4} = [];
    N{end,5} = mPeur > mPeuc;
    N{end,5} = N{end,5} + 1 - 1;
    N{end,6} = [];
    N{end,7} = [];
    N{end,8} = th_ue;
    N{end,9} = Wh_ue;
    N{end,10} = KSh_ue;
    N{end,11} = [];
    N{end,12} = mPuer > mPuec;
    N{end,12} = N{end,12} + 1 - 1;
    N{end,13} = [];
    N{end,14} = ClusterH;
    N{end,15} = ClusterM;
    Nname{end+1} = files{o}(1:end-8);
end

% Save
str = [{'t_hypothesis_eu'} {'W_hypothesis_eu'} {'KS_hypothesis_eu'}...
    {' '} {'eu_r>c'} {' '} {' '}...
    {'t_hypothesis_ue'} {'W_hypothesis_ue'} {'KS_hypothesis_ue'}...
    {' '} {'ue_r>c'} {' '} {'HC_cluster'} {'MS_cluster'}];
xlswrite(fn,str,'Noth','B1');
xlswrite(fn,Nname','Noth','A2');
xlswrite(fn,N,'Noth','B2');

save(fnm,'Noth','N')