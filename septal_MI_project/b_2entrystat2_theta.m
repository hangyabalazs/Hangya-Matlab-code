function [ThetaNoburst,ThetaBurst,TN,TB] = b_2entrystat2_theta
%2ENTRYSTAT2_THETA   Calculates statistics for entropy data (theta segments).
%   2ENTRYSTAT2_THETA runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   See also 2ENTRYSTAT2_NOTH and ENTRYSTAT_THETA.

% entropy
% theta
% eu r>c and ue r>c

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3_cont2\'];
inpplus = ['Theta\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3_cont2\Stat\power\theta_rorc.xls'];
fnm = [DATAPATH 'Entry3_cont2\Stat\power\stat_theta_rorc.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Theta\'];
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
ThetaNoburst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{},'t_hypothesis_ue',{},'W_hypothesis_ue',{},...
    'KS_hypothesis_ue',{});
ThetaBurst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{},'t_hypothesis_ue',{},'W_hypothesis_ue',{},...
    'KS_hypothesis_ue',{});
TN = {};
TB = {};
TNname = {};
TBname = {};
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
    load(files_clus{o})
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    ENT_eu_ctrl = aUxy;
    ENT_ue_ctrl = aUyx;
    mEeuc = mean(aUxy);
    mEuec = mean(aUyx);
    alfa = 0.007;
    [KSh_eu,KSp_eu] = kstest2(ENT_eu_real,ENT_eu_ctrl,alfa,'smaller');
    [KSh_ue,KSp_ue] = kstest2(ENT_ue_real,ENT_ue_ctrl,alfa,'smaller');
    [th_eu,tp_eu,ci_eu] = ttest(ENT_eu_real-ENT_eu_ctrl,0,alfa,'right');
    [th_ue,tp_ue,ci_ue] = ttest(ENT_ue_real-ENT_ue_ctrl,0,alfa,'right');
    [Wp_eu,Wh_eu] = b_signrank2(ENT_eu_real,ENT_eu_ctrl,'alpha',alfa);
    [Wp_ue,Wh_ue] = b_signrank2(ENT_ue_real,ENT_ue_ctrl,'alpha',alfa);
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
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-16);
        ThetaNoburst(end).t_hypothesis_eu = th_eu;
        ThetaNoburst(end).t_significance_eu = tp_eu;
        ThetaNoburst(end).W_hypothesis_eu = Wh_eu;
        ThetaNoburst(end).W_significance_eu = Wp_eu;
        ThetaNoburst(end).KS_hypothesis_eu = KSh_eu;
        ThetaNoburst(end).KS_significance_eu = KSp_eu;
        ThetaNoburst(end).t_hypothesis_ue = th_ue;
        ThetaNoburst(end).t_significance_ue = tp_ue;
        ThetaNoburst(end).W_hypothesis_ue = Wh_ue;
        ThetaNoburst(end).W_significance_ue = Wp_ue;
        ThetaNoburst(end).KS_hypothesis_ue = KSh_ue;
        ThetaNoburst(end).KS_significance_ue = KSp_ue;
        ThetaNoburst(end).mean_eu_real = mEeur;
        ThetaNoburst(end).mean_eu_ctrl = mEeuc;
        ThetaNoburst(end).mean_ue_real = mEuer;
        ThetaNoburst(end).mean_ue_ctrl = mEuec;
        TN{end+1,1} = th_eu;
        TN{end,2} = Wh_eu;
        TN{end,3} = KSh_eu;
        TN{end,4} = [];
        TN{end,5} = mEeur > mEeuc;
        TN{end,5} = TN{end,5} + 1 - 1;
        TN{end,6} = [];
        TN{end,7} = [];
        TN{end,8} = th_ue;
        TN{end,9} = Wh_ue;
        TN{end,10} = KSh_ue;
        TN{end,11} = [];
        TN{end,12} = mEuer > mEuec;
        TN{end,12} = TN{end,5} + 1 - 1;
        TNname{end+1} = files{o}(1:end-16);
    else
        ThetaBurst(end+1).name = files{o}(1:end-16);
        ThetaBurst(end).t_hypothesis_eu = th_eu;
        ThetaBurst(end).t_significance_eu = tp_eu;
        ThetaBurst(end).W_hypothesis_eu = Wh_eu;
        ThetaBurst(end).W_significance_eu = Wp_eu;
        ThetaBurst(end).KS_hypothesis_eu = KSh_eu;
        ThetaBurst(end).KS_significance_eu = KSp_eu;
        ThetaBurst(end).t_hypothesis_ue = th_ue;
        ThetaBurst(end).t_significance_ue = tp_ue;
        ThetaBurst(end).W_hypothesis_ue = Wh_ue;
        ThetaBurst(end).W_significance_ue = Wp_ue;
        ThetaBurst(end).KS_hypothesis_ue = KSh_ue;
        ThetaBurst(end).KS_significance_ue = KSp_ue;
        ThetaBurst(end).mean_eu_real = mEeur;
        ThetaBurst(end).mean_eu_ctrl = mEeuc;
        ThetaBurst(end).mean_ue_real = mEuer;
        ThetaBurst(end).mean_ue_ctrl = mEuec;
        TB{end+1,1} = th_eu;
        TB{end,2} = Wh_eu;
        TB{end,3} = KSh_eu;
        TB{end,4} = [];
        TB{end,5} = mEeur > mEeuc;
        TB{end,5} = TB{end,5} + 1 - 1;
        TB{end,6} = [];
        TB{end,7} = [];
        TB{end,8} = th_ue;
        TB{end,9} = Wh_ue;
        TB{end,10} = KSh_ue;
        TB{end,11} = [];
        TB{end,12} = mEuer > mEuec;
        TB{end,12} = TB{end,12} + 1 - 1;
        TBname{end+1} = files{o}(1:end-16);
    end
end
alfa

% Save
str = [{'t_hypothesis_eu'} {'W_hypothesis_eu'} {'KS_hypothesis_eu'}...
    {' '} {'eu_r>c'} {' '} {' '}...
    {'t_hypothesis_ue'} {'W_hypothesis_ue'} {'KS_hypothesis_ue'}...
    {' '} {'ue_r>c'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');

save(fnm,'ThetaNoburst','ThetaBurst','TN','TB')
cd(mm)