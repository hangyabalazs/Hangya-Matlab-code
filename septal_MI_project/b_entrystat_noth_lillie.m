function [NothNoburst,NothBurst,NN,NB] = b_entrystat_noth_lillie
%ENTRYSTAT_NOTH_LILLIE   Calculates statistics for entropy data (non-theta segments).
%   ENTRYSTAT_NOTH_LILLIE runs on the output of ENTRYRUN_THETA or ENTRYRUN3.
%   Edit code to specify exact direcctories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves Lilliefors-test
%   results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_THETA_LILLIEFORS.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry\'];
clusdir = [DATAPATH 'Burst\Cluster\Noth\'];
mm = pwd;
cd(inpdir)
cd(['Noth\dtf\halfsec\real'])
files0 = dir(pwd);
files = {};
files_clus = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if ext == '.mat'
        files{end+1} = files0(i).name;
        files_clus{end+1} = [files0(i).name(1:end-7) 'CLUSTER'];
    end
end
sf = length(files);

% Statistics
NothNoburst = struct('name',{},...
    'L_hypothesis_eu_real',{},'L_significance_eu_real',{},...
    'L_hypothesis_eu_ctrl',{},'L_significance_eu_ctrl',{},...
    'L_hypothesis_ue_real',{},'L_significance_ue_real',{},...
    'L_hypothesis_ue_ctrl',{},'L_significance_ue_ctrl',{});
NothBurst = struct('name',{},...
    'L_hypothesis_eu_real',{},'L_significance_eu_real',{},...
    'L_hypothesis_eu_ctrl',{},'L_significance_eu_ctrl',{},...
    'L_hypothesis_ue_real',{},'L_significance_ue_real',{},...
    'L_hypothesis_ue_ctrl',{},'L_significance_ue_ctrl',{});
NN = [];
NB = [];
NNname = {};
NBname = {};
for o = 1:sf
    cd(inpdir)
    cd(['Noth\dtf\halfsec\real'])
    load(files{o})
    PDC_eu_real = mean(PDC_eu(3:6,:));
    PDC_ue_real = mean(PDC_ue(3:6,:));
    mPeur = mean(PDC_eu_real);
    mPuer = mean(PDC_ue_real);
    mm2 = pwd;
    cd(clusdir)
    load(files_clus{o})
    cd(mm2)
    cd ..
    cd(['control'])
    load(files{o})
    PDC_eu_ctrl = mean(PDC_eu(3:6,:));
    PDC_ue_ctrl = mean(PDC_ue(3:6,:));
    mPeuc = mean(PDC_eu_ctrl);
    mPuec = mean(PDC_ue_ctrl);
    [Lh_eu_real,Lp_eu_real] = lillietest(PDC_eu_real);
    [Lh_eu_ctrl,Lp_eu_ctrl] = lillietest(PDC_eu_ctrl);
    [Lh_ue_real,Lp_ue_real] = lillietest(PDC_ue_real);
    [Lh_ue_ctrl,Lp_ue_ctrl] = lillietest(PDC_ue_ctrl);
    if ClusterNumber == 0
        NothNoburst(end+1).name = files{o}(1:end-8);
        NothNoburst(end).L_hypothesis_eu_real = Lh_eu_real;
        NothNoburst(end).L_significance_eu_real = Lp_eu_real;
        NothNoburst(end).L_hypothesis_eu_ctrl = Lh_eu_ctrl;
        NothNoburst(end).L_significance_eu_ctrl = Lp_eu_ctrl;
        NothNoburst(end).L_hypothesis_ue_real = Lh_ue_real;
        NothNoburst(end).L_significance_ue_real = Lp_ue_real;
        NothNoburst(end).L_hypothesis_ue_ctrl = Lh_ue_ctrl;
        NothNoburst(end).L_significance_ue_ctrl = Lp_ue_ctrl;
        NN(end+1,1) = Lh_eu_real;
        NN(end,2) = Lp_eu_real;
        NN(end,3) = Lh_eu_ctrl;
        NN(end,4) = Lp_eu_ctrl;
        NN(end,5) = Lh_ue_real;
        NN(end,6) = Lp_ue_real;
        NN(end,7) = Lh_ue_ctrl;
        NN(end,8) = Lp_ue_ctrl;
        NNname{end+1} = files{o}(1:end-8);
    else
        NothBurst(end+1).name = files{o}(1:end-8);
        NothBurst(end).L_hypothesis_eu_real = Lh_eu_real;
        NothBurst(end).L_significance_eu_real = Lp_eu_real;
        NothBurst(end).L_hypothesis_eu_ctrl = Lh_eu_ctrl;
        NothBurst(end).L_significance_eu_ctrl = Lp_eu_ctrl;
        NothBurst(end).L_hypothesis_ue_real = Lh_ue_real;
        NothBurst(end).L_significance_ue_real = Lp_ue_real;
        NothBurst(end).L_hypothesis_ue_ctrl = Lh_ue_ctrl;
        NothBurst(end).L_significance_ue_ctrl = Lp_ue_ctrl;
        NB(end+1,1) = Lh_eu_real;
        NB(end,2) = Lp_eu_real;
        NB(end,3) = Lh_eu_ctrl;
        NB(end,4) = Lp_eu_ctrl;
        NB(end,5) = Lh_ue_real;
        NB(end,6) = Lp_ue_real;
        NB(end,7) = Lh_ue_ctrl;
        NB(end,8) = Lp_ue_ctrl;
        NBname{end+1} = files{o}(1:end-8);
    end
end

% Save
fn = [DATAPATH 'Entry\Stat\noth_lillie.xls'];
str = [{'L_hypothesis_eu_real'} {'L_significance_eu_real'}...
    {'L_hypothesis_eu_ctrl'} {'L_significance_eu_ctrl'}...
    {'L_hypothesis_ue_real'} {'L_significance_ue_real'}...
    {'L_hypothesis_ue_ctrl'} {'L_significance_ue_ctrl'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');