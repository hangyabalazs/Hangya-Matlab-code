function [Noth,N] = b_entry3chstat4_noth
%ENTRY3CHSTAT4_NOTH   Calculates statistics for entropy data (theta segments).
%   ENTRY3CHSTAT4_NOTH runs on the output of ENTRYRUN_3CH. Edit code to
%   specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon ranksum-test (Mann-Whitney U-test) and
%   Kolmogorov-Smirnov test.
%
%   This program is for unit-unit analysis!
%
%   See also ENTRY3CHSTAT4_THETA and ENTRYSTAT4_NOTH.

% entropy
% noth
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry_unitunit\'];
inpplus = ['Noth\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry_unitunit\Stat\ENTROPY_power_ol1\2noth.xls'];
fnm = [DATAPATH 'Entry_unitunit\Stat\ENTROPY_power_ol1\stat_2noth.mat'];
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
Noth = struct([]);
N = {};
Nname = {};
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
    if mEeur > mEuer
        [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,0.05,'smaller');
        [th,tp,ci] = ttest2(ENT_eu_real,ENT_ue_real,[],'right','unequal');
        [Wp,Wh] = b_ranksum2(ENT_eu_real,ENT_ue_real);
    else
        [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,0.05,'larger');
        [th,tp,ci] = ttest2(ENT_eu_real,ENT_ue_real,[],'left','unequal');
        [Wp,Wh] = b_ranksum2(ENT_ue_real,ENT_eu_real);
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    Noth(end+1).name = files{o}(1:end-8);
    Noth(end).t_hypothesis = th;
    Noth(end).t_significance = tp;
    Noth(end).W_hypothesis = Wh;
    Noth(end).W_significance = Wp;
    Noth(end).KS_hypothesis = KSh;
    Noth(end).KS_significance = KSp;
    Noth(end).mean_eu_real = mEeur;
    Noth(end).mean_ue_real = mEuer;
    N{end+1,1} = th;
    N{end,2} = Wh;
    N{end,3} = KSh;
    N{end,4} = [];
    N{end,5} = mEeur > mEuer;
    N{end,5} = N{end,5} + 1 - 1;
    N{end,6} = [];
    N{end,7} = ClusterH;
    N{end,8} = ClusterM;
    Nname{end+1} = files{o}(1:end-8);
end

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'eu>?ue'} {' '} {'HC_cluster'} {'MS_cluster'}];
xlswrite(fn,str,'Noth','B1');
xlswrite(fn,Nname','Noth','A2');
xlswrite(fn,N,'Noth','B2');

save(fnm,'Noth','N')