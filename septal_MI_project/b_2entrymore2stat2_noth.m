function [NothNoburst,NothBurst,NN,NB] = b_2entrymore2stat2_noth
%2ENTRYMORE2STAT2_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRYMORE2STAT2_NOTH runs on the output of ENTRYRUN_NOTH or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   It analysis mutual information values instead of uncertainty coefficients.
%
%   See also 2ENTRYSTAT2_THETA and ENTRYSTAT_NOTH.

% entropy
% noth
% eu r>c and ue r>c

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3b_M\'];
inpplus = ['Noth\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3b_M\Stat\power\noth_mutinf_rorc.xls'];
fnm = [DATAPATH 'Entry3b_M\Stat\power\stat_noth_mutinf_rorc.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Noth\'];
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
NothNoburst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{});
NothBurst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{});
NN = {};
NB = {};
NNname = {};
NBname = {};
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
        NothNoburst(end+1).name = files{o}(1:end-16);
        NothNoburst(end).t_hypothesis_eu = th_eu;
        NothNoburst(end).t_significance_eu = tp_eu;
        NothNoburst(end).W_hypothesis_eu = Wh_eu;
        NothNoburst(end).W_significance_eu = Wp_eu;
        NothNoburst(end).KS_hypothesis_eu = KSh_eu;
        NothNoburst(end).KS_significance_eu = KSp_eu;
        NothNoburst(end).mean_eu_real = mEeur;
        NothNoburst(end).mean_eu_ctrl = mEeuc;
        NN{end+1,1} = th_eu;
        NN{end,2} = Wh_eu;
        NN{end,3} = KSh_eu;
        NN{end,4} = [];
        NN{end,5} = mEeur > mEeuc;
        NN{end,5} = NN{end,5} + 1 - 1;
        NNname{end+1} = files{o}(1:end-16);
    else
        NothBurst(end+1).name = files{o}(1:end-16);
        NothBurst(end).t_hypothesis_eu = th_eu;
        NothBurst(end).t_significance_eu = tp_eu;
        NothBurst(end).W_hypothesis_eu = Wh_eu;
        NothBurst(end).W_significance_eu = Wp_eu;
        NothBurst(end).KS_hypothesis_eu = KSh_eu;
        NothBurst(end).KS_significance_eu = KSp_eu;
        NothBurst(end).mean_eu_real = mEeur;
        NothBurst(end).mean_eu_ctrl = mEeuc;
        NB{end+1,1} = th_eu;
        NB{end,2} = Wh_eu;
        NB{end,3} = KSh_eu;
        NB{end,4} = [];
        NB{end,5} = mEeur > mEeuc;
        NB{end,5} = NB{end,5} + 1 - 1;;
        NBname{end+1} = files{o}(1:end-16);
    end
end
alfa

% Save
str = [{'t_hypothesis_Ieu'} {'W_hypothesis_Ieu'} {'KS_hypothesis_Ieu'}...
    {' '} {'Ieu_r>c'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');

save(fnm,'NothNoburst','NothBurst','NN','NB')
cd(mm)