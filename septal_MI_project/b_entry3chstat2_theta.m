function [Theta,T] = b_entry3chstat2_theta
%ENTRY3CHSTAT2_THETA   Calculates statistics for entropy data (theta segments).
%   ENTRY3CHSTAT2_THETA runs on the output of ENTRYRUN_3CH. Edit code to
%   specify exact directories!
%
%   It runs on entropy results, and addresses the question whether real data is
%   larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon ranksum-test (Mann-Whitney U-test) and Kolmogorov-
%   Smirnov test.
%
%   This program is for unit-unit analysis!
%
%   See also ENTRY3CHSTAT2_NOTH and ENTRYSTAT2_THETA.

% entropy
% theta
% eu r>c and ue r>c

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry2_unitunit\'];
inpplus = ['Theta\entropy\phase\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry2_unitunit\Stat\ENTROPY_phase_ol1\theta_rorc.xls'];
fnm = [DATAPATH 'Entry2_unitunit\Stat\ENTROPY_phase_ol1\stat_theta_rorc.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Theta_3ch\'];
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
        files_clus{end+1} = [files0(i).name(1:end-15) 'CLUSTER'];
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
    ENT_eu_real = aUxy;
    ENT_ue_real = aUyx;
    mEeur = mean(aUxy);
    mEuer = mean(aUyx);
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
    ENT_eu_ctrl = aUxy;
    ENT_ue_ctrl = aUyx;
    mEeuc = mean(aUxy);
    mEuec = mean(aUyx);
    [KSh_eu,KSp_eu] = kstest2(ENT_eu_real,ENT_eu_ctrl,0.05,'smaller');
    [KSh_ue,KSp_ue] = kstest2(ENT_ue_real,ENT_ue_ctrl,0.05,'smaller');
    [th_eu,tp_eu,ci_eu] = ttest2(ENT_eu_real,ENT_eu_ctrl,[],'right','unequal');
    [th_ue,tp_ue,ci_ue] = ttest2(ENT_ue_real,ENT_ue_ctrl,[],'right','unequal');
    [Wp_eu,Wh_eu] = b_ranksum2(ENT_eu_real,ENT_eu_ctrl);
    [Wp_ue,Wh_ue] = b_ranksum2(ENT_ue_real,ENT_ue_ctrl);
    KSh_eu = KSh_eu + 1 - 1;    % convert to numbers from logical values
    th_eu = th_eu + 1 - 1;
    Wh_eu = Wh_eu + 1 - 1;
    KSh_ue = KSh_ue + 1 - 1;
    th_ue = th_ue + 1 - 1;
    Wh_ue = Wh_ue + 1 - 1;
    if mEeuc > mEeur
        Wh_eu = 0;
        Wp_eu = -1;
    end
    if mEuec > mEuer
        Wh_ue = 0;
        Wp_ue = -1;
    end
    Theta(end+1).name = files{o}(1:end-16);
    Theta(end).t_hypothesis_eu = th_eu;
    Theta(end).t_significance_eu = tp_eu;
    Theta(end).W_hypothesis_eu = Wh_eu;
    Theta(end).W_significance_eu = Wp_eu;
    Theta(end).KS_hypothesis_eu = KSh_eu;
    Theta(end).KS_significance_eu = KSp_eu;
    Theta(end).t_hypothesis_ue = th_ue;
    Theta(end).t_significance_ue = tp_ue;
    Theta(end).W_hypothesis_ue = Wh_ue;
    Theta(end).W_significance_ue = Wp_ue;
    Theta(end).KS_hypothesis_ue = KSh_ue;
    Theta(end).KS_significance_ue = KSp_ue;
    Theta(end).mean_eu_real = mEeur;
    Theta(end).mean_eu_ctrl = mEeuc;
    Theta(end).mean_ue_real = mEuer;
    Theta(end).mean_ue_ctrl = mEuec;
    T{end+1,1} = th_eu;
    T{end,2} = Wh_eu;
    T{end,3} = KSh_eu;
    T{end,4} = [];
    T{end,5} = mEeur > mEeuc;
    T{end,5} = T{end,5} + 1 - 1;
    T{end,6} = [];
    T{end,7} = [];
    T{end,8} = th_ue;
    T{end,9} = Wh_ue;
    T{end,10} = KSh_ue;
    T{end,11} = [];
    T{end,12} = mEuer > mEuec;
    T{end,12} = T{end,5} + 1 - 1;
    T{end,13} = [];
    T{end,14} = ClusterH;
    T{end,15} = ClusterM;
    Tname{end+1} = files{o}(1:end-16);
end

% Save
str = [{'t_hypothesis_eu'} {'W_hypothesis_eu'} {'KS_hypothesis_eu'}...
    {' '} {'eu_r>c'} {' '} {' '}...
    {'t_hypothesis_ue'} {'W_hypothesis_ue'} {'KS_hypothesis_ue'}...
    {' '} {'ue_r>c'} {' '} {'HC_cluster'} {'MS_cluster'}];
xlswrite(fn,str,'Theta','B1');
xlswrite(fn,Tname','Theta','A2');
xlswrite(fn,T,'Theta','B2');

save(fnm,'Theta','T')