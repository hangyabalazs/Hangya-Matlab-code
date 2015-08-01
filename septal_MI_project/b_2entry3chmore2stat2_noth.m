function [Noth,N] = b_2entry3chmore2stat2_noth
%2ENTRY3CHMORE2STAT2_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRY3CHMORE2STAT2_NOTH runs on the output of ENTRYRUN_NOTH or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   It analysis mutual information values instead of uncertainty coefficients.
%
%   This program is for unit-unit analysis!
%
%   See also 2ENTRYSTAT2_THETA and ENTRYSTAT_NOTH.

% entropy
% noth
% eu r>c and ue r>c

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3c_unitunit\'];
inpplus = ['Noth\entropy\phase\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3c_unitunit\Stat\phase\noth_mutinf_rorc.xls'];
fnm = [DATAPATH 'Entry3c_unitunit\Stat\phase\stat_noth_mutinf_rorc.mat'];
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
Noth = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{});
N = {};
Nname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    ENT_eu_real = aIxy;
    mEeur = mean(aIxy);
    mm2 = pwd;
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
    Noth(end+1).name = files{o}(1:end-16);
    Noth(end).t_hypothesis_eu = th_eu;
    Noth(end).t_significance_eu = tp_eu;
    Noth(end).W_hypothesis_eu = Wh_eu;
    Noth(end).W_significance_eu = Wp_eu;
    Noth(end).KS_hypothesis_eu = KSh_eu;
    Noth(end).KS_significance_eu = KSp_eu;
    Noth(end).mean_eu_real = mEeur;
    Noth(end).mean_eu_ctrl = mEeuc;
    N{end+1,1} = th_eu;
    N{end,2} = Wh_eu;
    N{end,3} = KSh_eu;
    N{end,4} = [];
    N{end,5} = mEeur > mEeuc;
    N{end,5} = N{end,5} + 1 - 1;
    Nname{end+1} = files{o}(1:end-16);
end
alfa

% Save
str = [{'t_hypothesis_Ieu'} {'W_hypothesis_Ieu'} {'KS_hypothesis_Ieu'}...
    {' '} {'Ieu_r>c'}];
xlswrite(fn,str,'Noth','B1');
xlswrite(fn,Nname','Noth','A2');
xlswrite(fn,N,'Noth','B2');

save(fnm,'Noth','N')
cd(mm)