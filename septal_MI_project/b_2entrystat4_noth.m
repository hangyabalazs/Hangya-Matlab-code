function [NothNoburst,NothBurst,NN,NB] = b_2entrystat4_noth
%2ENTRYSTAT4_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRYSTAT4_NOTH runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   See also 2ENTRYSTAT4_THETA and ENTRYSTAT3_NOTH.

% entropy
% noth
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3_cont2\'];
inpplus = ['Noth\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3_cont2\Stat\power\2noth.xls'];
fnm = [DATAPATH 'Entry3_cont2\Stat\power\stat_2noth.mat'];
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
NothNoburst = struct('name',{},'t_hypothesis',{},'W_hypothesis',{},...
    'KS_hypothesis',{});
NothBurst = struct('name',{},'t_hypothesis',{},'W_hypothesis',{},...
    'KS_hypothesis',{});
NN = {};
NB = {};
NNname = {};
NBname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
    load(files{o})
    aUxy = round(aUxy*1000) / 1000;
    aUyx = round(aUyx*1000) / 1000;
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
    if ClusterNumber == 0
        NothNoburst(end+1).name = files{o}(1:end-8);
        NothNoburst(end).t_hypothesis = th;
        NothNoburst(end).t_significance = tp;
        NothNoburst(end).W_hypothesis = Wh;
        NothNoburst(end).W_significance = Wp;
        NothNoburst(end).KS_hypothesis = KSh;
        NothNoburst(end).KS_significance = KSp;
        NothNoburst(end).mean_eu_real = mEeur;
        NothNoburst(end).mean_ue_real = mEuer;
        NN{end+1,1} = th;
        NN{end,2} = Wh;
        NN{end,3} = KSh;
        NN{end,4} = [];
        NN{end,5} = mEeur > mEuer;
        NN{end,5} = NN{end,5} + 1 - 1;
        NNname{end+1} = files{o}(1:end-8);
    else
        NothBurst(end+1).name = files{o}(1:end-8);
        NothBurst(end).t_hypothesis = th;
        NothBurst(end).t_significance = tp;
        NothBurst(end).W_hypothesis = Wh;
        NothBurst(end).W_significance = Wp;
        NothBurst(end).KS_hypothesis = KSh;
        NothBurst(end).KS_significance = KSp;
        NothBurst(end).mean_eu_real = mEeur;
        NothBurst(end).mean_ue_real = mEuer;
        NB{end+1,1} = th;
        NB{end,2} = Wh;
        NB{end,3} = KSh;
        NB{end,4} = [];
        NB{end,5} = mEeur > mEuer;
        NB{end,5} = NB{end,5} + 1 - 1;
        NBname{end+1} = files{o}(1:end-8);
    end
end
alfa

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'eu>?ue'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');

save(fnm,'NothNoburst','NothBurst','NN','NB')
cd(mm)