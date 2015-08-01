function [Theta,T] = b_2entry3chstat4_theta
%2ENTRY3CHSTAT4_THETA   Calculates statistics for entropy data (theta segments).
%   2ENTRY3CHSTAT4_THETA runs on the output of ENTRYRUN_3CH. Edit code to
%   specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   This program is for unit-unit analysis!
%
%   See also 2ENTRY3CHSTAT4_NOTH and ENTRYSTAT4_THETA.

% entropy
% theta
% eu>?ue

%Directories
name = '';
global DATAPATH
inpdir = [DATAPATH 'Entry3c_cont3_unitunit\'];
inpplus = ['Theta\entropy\phase\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3c_cont3_unitunit\Stat\phase\2theta.xls'];
fnm = [DATAPATH 'Entry3c_cont3_unitunit\Stat\phase\stat_2theta.mat'];
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
    alfa = 0.007;
    if mEeur > mEuer
        [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,alfa,'smaller');
        [th,tp,ci] = ttest(ENT_eu_real-ENT_ue_real,0,alfa,'right');
        [Wp,Wh] = b_signrank2(ENT_eu_real,ENT_ue_real,'alpha',alfa);
    else
        [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,alfa,'larger');
        [th,tp,ci] = ttest(ENT_eu_real-ENT_ue_real,0,alfa,'left');
        [Wp,Wh] = b_signrank2(ENT_ue_real,ENT_eu_real,'alpha',alfa);
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    Theta(end+1).name = files{o}(1:end-8);
    Theta(end).t_hypothesis = th;
    Theta(end).t_significance = tp;
    Theta(end).W_hypothesis = Wh;
    Theta(end).W_significance = Wp;
    Theta(end).KS_hypothesis = KSh;
    Theta(end).KS_significance = KSp;
    Theta(end).mean_eu_real = mEeur;
    Theta(end).mean_ue_real = mEuer;
    Theta(end).HC_unit_cluternumber = ClusterH;
    Theta(end).MS_unit_cluternumber = ClusterM;
    T{end+1,1} = th;
    T{end,2} = Wh;
    T{end,3} = KSh;
    T{end,4} = [];
    T{end,5} = mEeur > mEuer;
    T{end,5} = T{end,5} + 1 - 1;
    T{end,6} = [];
    T{end,7} = ClusterH;
    T{end,8} = ClusterM;
    Tname{end+1} = files{o}(1:end-8);
end
alfa

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'eu>?ue'} {' '} {'HC_cluster'} {'MS_cluster'}];
xlswrite(fn,str,'Theta','B1');
xlswrite(fn,Tname','Theta','A2');
xlswrite(fn,T,'Theta','B2');

save(fnm,'Theta','T')
cd(mm)