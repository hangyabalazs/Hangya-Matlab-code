function [Theta,T] = b_2entry3chmmstat4_theta
%2ENTRY3CHMMSTAT4_THETA   Calculates statistics for entropy data (theta segments).
%   2ENTRY3CHMMSTAT4_THETA runs on the output of ENTRYRUN_THETA or ENTRYRUN3C. Edit
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
%   See also 2ENTRYSTAT4_NOTH and ENTRYSTAT3_THETA.

% entropy
% theta
% eu>?ue

%Directories
name = '';
global DATAPATH
inpdir = [DATAPATH 'Entry3c_unitunit\'];
inpplus = ['Theta\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3c_unitunit\Stat\power\2theta_eu.xls'];
fnm = [DATAPATH 'Entry3c_unitunit\Stat\power\stat_2theta_eu.mat'];
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
Theta = struct('name',{},'t_hypothesis_diff',{},'W_hypothesis_diff',{},...
    'KS_hypothesis_diff',{});
T = {};
Tname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    ENT_eu_real = aUxy;
    ENT_ue_real = aUyx;
    ENT_diff_real = ENT_eu_real - ENT_ue_real;
    mEdiffr = mean(ENT_diff_real);
    mEuer = mean(aUyx);
    mm2 = pwd;
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    ENT_eu_ctrl = aUxy;
    ENT_ue_ctrl = aUyx;
    ENT_diff_ctrl = ENT_eu_ctrl - ENT_ue_ctrl;
    mEdiffc = mean(ENT_diff_ctrl);
    mEuec = mean(aUyx);
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
    Theta(end+1).name = files{o}(1:end-16);
    Theta(end).t_hypothesis_diff = th_diff;
    Theta(end).t_significance_diff = tp_diff;
    Theta(end).W_hypothesis_diff = Wh_diff;
    Theta(end).W_significance_diff = Wp_diff;
    Theta(end).KS_hypothesis_diff = KSh_diff;
    Theta(end).KS_significance_diff = KSp_diff;
    Theta(end).mean_diff_real = mEdiffr;
    Theta(end).mean_diff_ctrl = mEdiffc;
    T{end+1,1} = th_diff;
    T{end,2} = Wh_diff;
    T{end,3} = KSh_diff;
    T{end,4} = [];
    T{end,5} = mEdiffr > mEdiffc;
    T{end,5} = T{end,5} + 1 - 1;
    Tname{end+1} = files{o}(1:end-16);
end
alfa

% Save
str = [{'t_hypothesis_diff'} {'W_hypothesis_diff'} {'KS_hypothesis_diff'}...
    {' '} {'diff_r>c'}];
xlswrite(fn,str,'Theta','B1');
xlswrite(fn,Tname','Theta','A2');
xlswrite(fn,T,'Theta','B2');

save(fnm,'Theta','T')
cd(mm)