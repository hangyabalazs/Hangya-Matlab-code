function [ThetaNoburst,ThetaBurst,TN,TB] = b_entrystat_theta_lillie
%ENTRYSTAT_THETA_LILLIE   Calculates statistics for entropy data (theta segments).
%   ENTRYSTAT_THETA_LILLIE runs on the output of ENTRYRUN_THETA or ENTRYRUN3.
%   Edit code to specify exact direcctories!
%
%   It runs on PDC results, and addresses the question whether real data
%   is larger than control or not. It calculates and saves Lilliefors-test
%   results.
%
%   Not used for final analysis!
%
%   See also ENTRYSTAT_NOTH_LILLIEFORS.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry\'];
clusdir = [DATAPATH 'Burst\Cluster\Theta\'];
mm = pwd;
cd(inpdir)
cd(['Theta\dtf\halfsec\real'])
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
ThetaNoburst = struct('name',{},...
    'L_hypothesis_eu_real',{},'L_significance_eu_real',{},...
    'L_hypothesis_eu_ctrl',{},'L_significance_eu_ctrl',{},...
    'L_hypothesis_ue_real',{},'L_significance_ue_real',{},...
    'L_hypothesis_ue_ctrl',{},'L_significance_ue_ctrl',{});
ThetaBurst = struct('name',{},...
    'L_hypothesis_eu_real',{},'L_significance_eu_real',{},...
    'L_hypothesis_eu_ctrl',{},'L_significance_eu_ctrl',{},...
    'L_hypothesis_ue_real',{},'L_significance_ue_real',{},...
    'L_hypothesis_ue_ctrl',{},'L_significance_ue_ctrl',{});
TN = [];
TB = [];
TNname = {};
TBname = {};
for o = 1:sf
    cd(inpdir)
    cd(['Theta\dtf\halfsec\real'])
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
        ThetaNoburst(end+1).name = files{o}(1:end-8);
        ThetaNoburst(end).L_hypothesis_eu_real = Lh_eu_real;
        ThetaNoburst(end).L_significance_eu_real = Lp_eu_real;
        ThetaNoburst(end).L_hypothesis_eu_ctrl = Lh_eu_ctrl;
        ThetaNoburst(end).L_significance_eu_ctrl = Lp_eu_ctrl;
        ThetaNoburst(end).L_hypothesis_ue_real = Lh_ue_real;
        ThetaNoburst(end).L_significance_ue_real = Lp_ue_real;
        ThetaNoburst(end).L_hypothesis_ue_ctrl = Lh_ue_ctrl;
        ThetaNoburst(end).L_significance_ue_ctrl = Lp_ue_ctrl;
        TN(end+1,1) = Lh_eu_real;
        TN(end,2) = Lp_eu_real;
        TN(end,3) = Lh_eu_ctrl;
        TN(end,4) = Lp_eu_ctrl;
        TN(end,5) = Lh_ue_real;
        TN(end,6) = Lp_ue_real;
        TN(end,7) = Lh_ue_ctrl;
        TN(end,8) = Lp_ue_ctrl;
        TNname{end+1} = files{o}(1:end-8);
    else
        ThetaBurst(end+1).name = files{o}(1:end-8);
        ThetaBurst(end).L_hypothesis_eu_real = Lh_eu_real;
        ThetaBurst(end).L_significance_eu_real = Lp_eu_real;
        ThetaBurst(end).L_hypothesis_eu_ctrl = Lh_eu_ctrl;
        ThetaBurst(end).L_significance_eu_ctrl = Lp_eu_ctrl;
        ThetaBurst(end).L_hypothesis_ue_real = Lh_ue_real;
        ThetaBurst(end).L_significance_ue_real = Lp_ue_real;
        ThetaBurst(end).L_hypothesis_ue_ctrl = Lh_ue_ctrl;
        ThetaBurst(end).L_significance_ue_ctrl = Lp_ue_ctrl;
        TB(end+1,1) = Lh_eu_real;
        TB(end,2) = Lp_eu_real;
        TB(end,3) = Lh_eu_ctrl;
        TB(end,4) = Lp_eu_ctrl;
        TB(end,5) = Lh_ue_real;
        TB(end,6) = Lp_ue_real;
        TB(end,7) = Lh_ue_ctrl;
        TB(end,8) = Lp_ue_ctrl;
        TBname{end+1} = files{o}(1:end-8);
    end
end

% Save
fn = [DATAPATH 'Entry\Stat\theta_lillie.xls'];
str = [{'L_hypothesis_eu_real'} {'L_significance_eu_real'}...
    {'L_hypothesis_eu_ctrl'} {'L_significance_eu_ctrl'}...
    {'L_hypothesis_ue_real'} {'L_significance_ue_real'}...
    {'L_hypothesis_ue_ctrl'} {'L_significance_ue_ctrl'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');