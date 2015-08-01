function [ThetaNoburst,ThetaBurst,TN,TB] = b_2entrymore2stat2_theta
%2ENTRYMORE2STAT2_THETA   Calculates statistics for entropy data (theta segments).
%   2ENTRYMORE2STAT2_THETA runs on the output of ENTRYRUN_THETA or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   It analysis mutual information values instead of uncertainty coefficients.
%
%   See also 2ENTRYSTAT2_NOTH and ENTRYSTAT_THETA.

% entropy
% theta
% eu r>c and ue r>c

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3b_M\'];
inpplus = ['Theta\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3b_M\Stat\power\theta_mutinf_rorc.xls'];
fnm = [DATAPATH 'Entry3b_M\Stat\power\stat_theta_mutinf_rorc.mat'];
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
    'KS_hypothesis_eu',{});
ThetaBurst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{});
TN = {};
TB = {};
TNname = {};
TBname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    ENT_eu_real = aIxy;
    mEeur = mean(aIxy);
    mm2 = pwd;
    cd(clusdir)
    load(files_clus{o})
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    ENT_eu_ctrl = aIxy;
    mEeuc = mean(aIxy);
    alfa = 0.007;
    [KSh_eu,KSp_eu] = kstest2(ENT_eu_real,ENT_eu_ctrl,alfa,'smaller');
    [th_eu,tp_eu,ci_eu] = ttest(ENT_eu_real-ENT_eu_ctrl,0,alfa,'right');
    [Wp_eu,Wh_eu] = b_signrank2(ENT_eu_real,ENT_eu_ctrl,'alpha',alfa);
    KSh_eu = KSh_eu + 1 - 1;    % convert to numbers from logical values
    th_eu = th_eu + 1 - 1;
    Wh_eu = Wh_eu + 1 - 1;
    if mEeuc > mEeur
        Wh_eu = 0;
        Wp_eu = -1;
    end
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-16);
        ThetaNoburst(end).t_hypothesis_eu = th_eu;
        ThetaNoburst(end).t_significance_eu = tp_eu;
        ThetaNoburst(end).W_hypothesis_eu = Wh_eu;
        ThetaNoburst(end).W_significance_eu = Wp_eu;
        ThetaNoburst(end).KS_hypothesis_eu = KSh_eu;
        ThetaNoburst(end).KS_significance_eu = KSp_eu;
        ThetaNoburst(end).mean_eu_real = mEeur;
        ThetaNoburst(end).mean_eu_ctrl = mEeuc;
        TN{end+1,1} = th_eu;
        TN{end,2} = Wh_eu;
        TN{end,3} = KSh_eu;
        TN{end,4} = [];
        TN{end,5} = mEeur > mEeuc;
        TN{end,5} = TN{end,5} + 1 - 1;
        TNname{end+1} = files{o}(1:end-16);
    else
        ThetaBurst(end+1).name = files{o}(1:end-16);
        ThetaBurst(end).t_hypothesis_eu = th_eu;
        ThetaBurst(end).t_significance_eu = tp_eu;
        ThetaBurst(end).W_hypothesis_eu = Wh_eu;
        ThetaBurst(end).W_significance_eu = Wp_eu;
        ThetaBurst(end).KS_hypothesis_eu = KSh_eu;
        ThetaBurst(end).KS_significance_eu = KSp_eu;
        ThetaBurst(end).mean_eu_real = mEeur;
        ThetaBurst(end).mean_eu_ctrl = mEeuc;
        TB{end+1,1} = th_eu;
        TB{end,2} = Wh_eu;
        TB{end,3} = KSh_eu;
        TB{end,4} = [];
        TB{end,5} = mEeur > mEeuc;
        TB{end,5} = TB{end,5} + 1 - 1;
        TBname{end+1} = files{o}(1:end-16);
    end
end
alfa

% Save
str = [{'t_hypothesis_Ieu'} {'W_hypothesis_Ieu'} {'KS_hypothesis_Ieu'}...
    {' '} {'Ieu_r>c'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');

save(fnm,'ThetaNoburst','ThetaBurst','TN','TB')
cd(mm)