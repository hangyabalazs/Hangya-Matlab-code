function [NothNoburst,NothBurst,NN,NB] = b_2entrystat_noth
%2ENTRYSTAT_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRYSTAT_NOTH runs on the output of ENTRYRUN_NOTH or ENTRYRUN3. Edit
%   code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether real data is
%   larger than control or not. It calculates and saves 3 statistics:
%   t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   See also 2ENTRYSTAT_THETA and ENTRYSTAT2_NOTH.

% DTF
% noth
% eu r>c and ue r>c

% Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3_cont2\'];
inpplus = ['Noth\dtf\onesec\real'];
fn = [DATAPATH 'Entry3_cont2\Stat\DTF_onesec\noth_rorc.xls'];
fnm = [DATAPATH 'Entry3_cont2\Stat\DTF_onesec\stat_noth_rorc.mat'];
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
        files_clus{end+1} = [files0(i).name(1:end-7) 'CLUSTER'];
    end
end
sf = length(files);

% Statistics
NothNoburst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{},'t_hypothesis_ue',{},'W_hypothesis_ue',{},...
    'KS_hypothesis_ue',{});
NothBurst = struct('name',{},'t_hypothesis_eu',{},'W_hypothesis_eu',{},...
    'KS_hypothesis_eu',{},'t_hypothesis_ue',{},'W_hypothesis_ue',{},...
    'KS_hypothesis_ue',{});
NN = {};
NB = {};
NNname = {};
NBname = {};
for o = 1:sf
    cd(inpdir)
    cd(inpplus)
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
    alfa = 0.007;
    [KSh_eu,KSp_eu] = kstest2(PDC_eu_real,PDC_eu_ctrl,alfa,'smaller');
    [KSh_ue,KSp_ue] = kstest2(PDC_ue_real,PDC_ue_ctrl,alfa,'smaller');
    [th_eu,tp_eu,ci_eu] = ttest(PDC_eu_real-PDC_eu_ctrl,0,alfa,'right');
    [th_ue,tp_ue,ci_ue] = ttest(PDC_ue_real-PDC_ue_ctrl,0,alfa,'right');
    [Wp_eu,Wh_eu] = b_signrank2(PDC_eu_real,PDC_eu_ctrl,'alpha',alfa);
    [Wp_ue,Wh_ue] = b_signrank2(PDC_ue_real,PDC_ue_ctrl,'alpha',alfa);
    KSh_eu = KSh_eu + 1 - 1;    % convert to numbers from logical values
    th_eu = th_eu + 1 - 1;
    Wh_eu = Wh_eu + 1 - 1;
    KSh_ue = KSh_ue + 1 - 1;
    th_ue = th_ue + 1 - 1;
    Wh_ue = Wh_ue + 1 - 1;
    if mPeuc > mPeur
        Wh_eu = 0;
        Wp_eu = -1;
    end
    if mPuec > mPuer
        Wh_ue = 0;
        Wp_ue = -1;
    end
    if ClusterNumber == 0
        NothNoburst(end+1).name = files{o}(1:end-8);
        NothNoburst(end).t_hypothesis_eu = th_eu;
        NothNoburst(end).t_significance_eu = tp_eu;
        NothNoburst(end).W_hypothesis_eu = Wh_eu;
        NothNoburst(end).W_significance_eu = Wp_eu;
        NothNoburst(end).KS_hypothesis_eu = KSh_eu;
        NothNoburst(end).KS_significance_eu = KSp_eu;
        NothNoburst(end).t_hypothesis_ue = th_ue;
        NothNoburst(end).t_significance_ue = tp_ue;
        NothNoburst(end).W_hypothesis_ue = Wh_ue;
        NothNoburst(end).W_significance_ue = Wp_ue;
        NothNoburst(end).KS_hypothesis_ue = KSh_ue;
        NothNoburst(end).KS_significance_ue = KSp_ue;
        NothNoburst(end).mean_eu_real = mPeur;
        NothNoburst(end).mean_eu_ctrl = mPeuc;
        NothNoburst(end).mean_ue_real = mPuer;
        NothNoburst(end).mean_ue_ctrl = mPuec;
        NN{end+1,1} = th_eu;
        NN{end,2} = Wh_eu;
        NN{end,3} = KSh_eu;
        NN{end,4} = [];
        NN{end,5} = mPeur > mPeuc;
        NN{end,5} = NN{end,5} + 1 - 1;
        NN{end,6} = [];
        NN{end,7} = [];
        NN{end,8} = th_ue;
        NN{end,9} = Wh_ue;
        NN{end,10} = KSh_ue;
        NN{end,11} = [];
        NN{end,12} = mPuer > mPuec;
        NN{end,12} = NN{end,12} + 1 - 1;
        NNname{end+1} = files{o}(1:end-8);
    else
        NothBurst(end+1).name = files{o}(1:end-8);
        NothBurst(end).t_hypothesis_eu = th_eu;
        NothBurst(end).t_significance_eu = tp_eu;
        NothBurst(end).W_hypothesis_eu = Wh_eu;
        NothBurst(end).W_significance_eu = Wp_eu;
        NothBurst(end).KS_hypothesis_eu = KSh_eu;
        NothBurst(end).KS_significance_eu = KSp_eu;
        NothBurst(end).t_hypothesis_ue = th_ue;
        NothBurst(end).t_significance_ue = tp_ue;
        NothBurst(end).W_hypothesis_ue = Wh_ue;
        NothBurst(end).W_significance_ue = Wp_ue;
        NothBurst(end).KS_hypothesis_ue = KSh_ue;
        NothBurst(end).KS_significance_ue = KSp_ue;
        NothBurst(end).mean_eu_real = mPeur;
        NothBurst(end).mean_eu_ctrl = mPeuc;
        NothBurst(end).mean_ue_real = mPuer;
        NothBurst(end).mean_ue_ctrl = mPuec;
        NB{end+1,1} = th_eu;
        NB{end,2} = Wh_eu;
        NB{end,3} = KSh_eu;
        NB{end,4} = [];
        NB{end,5} = mPeur > mPeuc;
        NB{end,5} = NB{end,5} + 1 - 1;;
        NB{end,6} = [];
        NB{end,7} = [];
        NB{end,8} = th_ue;
        NB{end,9} = Wh_ue;
        NB{end,10} = KSh_ue;
        NB{end,11} = [];
        NB{end,12} = mPuer > mPuec;
        NB{end,12} = NB{end,12} + 1 - 1;
        NBname{end+1} = files{o}(1:end-8);
    end
end
alfa

% Save
str = [{'t_hypothesis_eu'} {'W_hypothesis_eu'} {'KS_hypothesis_eu'}...
    {' '} {'eu_r>c'} {' '} {' '}...
    {'t_hypothesis_ue'} {'W_hypothesis_ue'} {'KS_hypothesis_ue'}...
    {' '} {'ue_r>c'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');

save(fnm,'NothNoburst','NothBurst','NN','NB')
cd(mm)