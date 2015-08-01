function [Noth,N] = b_2entry3chmmstat4_noth
%2ENTRY3CHMMSTAT4_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRY3CHMMSTAT4_NOTH runs on the output of ENTRYRUN_THETA or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS* 
%   information flow is larger than control HC-MS or not. It calculates and
%   saves 3 statistics: t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov
%   test.
%   *: "eu" in the result file name indicates HC-MS, whereas "ue" stands
%   for MS-HC. To convert from "eu" to "ue",
%       ENT_diff_real = ENT_eu_real - ENT_ue_real;
%       ENT_diff_ctrl = ENT_eu_ctrl - ENT_ue_ctrl;
%   rows should be modified to
%       ENT_diff_real = ENT_ue_real - ENT_eu_real;
%       ENT_diff_ctrl = ENT_ue_ctrl - ENT_eu_ctrl;
%
%   See also 2ENTRYSTAT4_THETA and ENTRYSTAT3_NOTH.

% entropy
% noth
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3c_unitunit\'];
inpplus = ['Noth\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3c_unitunit\Stat\power\2noth_eu.xls'];
fnm = [DATAPATH 'Entry3c_unitunit\Stat\power\stat_2noth_eu.mat'];
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
        files_clus{end+1} = [files0(i).name(1:end-15) 'CLUSTER'];
    end
end
sf = length(files);

% Statistics
Noth = struct('name',{},'t_hypothesis_diff',{},'W_hypothesis_diff',{},...
    'KS_hypothesis_diff',{});
N = {};
Nname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    ENT_eu_real = aUxy;
    ENT_ue_real = aUyx;
    ENT_diff_real = ENT_eu_real - ENT_ue_real;
    mEdiffr = mean(ENT_diff_real);
    mm2 = pwd;
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    ENT_eu_ctrl = aUxy;
    ENT_ue_ctrl = aUyx;
    ENT_diff_ctrl = ENT_eu_ctrl - ENT_ue_ctrl;
    mEdiffc = mean(ENT_diff_ctrl);
    alfa = 0.007;
    [KSh_diff,KSp_diff] = kstest2(ENT_diff_real,ENT_diff_ctrl,alfa,'smaller');
    [th_diff,tp_diff,ci_diff] = ttest(ENT_diff_real-ENT_diff_ctrl,0,alfa,'right');
    [Wp_diff,Wh_diff] = b_signrank2(ENT_diff_real,ENT_diff_ctrl,'alpha',alfa);
    KSh_diff = KSh_diff + 1 - 1;    % convert to numbers from logical values
    th_diff = th_diff + 1 - 1;
    Wh_diff = Wh_diff + 1 - 1;
    if mEdiffc > mEdiffr
        Wh_diff = 0;
        Wp_diff = -1;
    end
    Noth(end+1).name = files{o}(1:end-16);
    Noth(end).t_hypothesis_diff = th_diff;
    Noth(end).t_significance_diff = tp_diff;
    Noth(end).W_hypothesis_diff = Wh_diff;
    Noth(end).W_significance_diff = Wp_diff;
    Noth(end).KS_hypothesis_diff = KSh_diff;
    Noth(end).KS_significance_diff = KSp_diff;
    Noth(end).mean_diff_real = mEdiffr;
    Noth(end).mean_diff_ctrl = mEdiffc;
    N{end+1,1} = th_diff;
    N{end,2} = Wh_diff;
    N{end,3} = KSh_diff;
    N{end,4} = [];
    N{end,5} = mEdiffr > mEdiffc;
    N{end,5} = N{end,5} + 1 - 1;
    Nname{end+1} = files{o}(1:end-16);
end
alfa

% Save
str = [{'t_hypothesis_diff'} {'W_hypothesis_diff'} {'KS_hypothesis_diff'}...
    {' '} {'diff_r>c'}];
xlswrite(fn,str,'Noth','B1');
xlswrite(fn,Nname','Noth','A2');
xlswrite(fn,N,'Noth','B2');

save(fnm,'Noth','N')
cd(mm)