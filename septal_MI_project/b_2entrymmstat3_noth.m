function [NothNoburst,NothBurst,NN,NB] = b_2entrymmstat3_noth
%2ENTRYMMSTAT3_NOTH   Calculates statistics for entropy data (non-theta segments).
%   2ENTRYMMSTAT3_NOTH runs on the output of ENTRYRUN_THETA or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on PDC results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon signrank-test and Kolmogorov-Smirnov test.
%
%   It compares real ue-eu difference to control difference.
%
%   See also 2ENTRYSTAT3_THETA and ENTRYSTAT4_NOTH.

% DTF
% noth
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry3b\'];
inpplus = ['Noth\dtf\onesec\real\'];
fn = [DATAPATH 'Entry3b_GG\Stat\DTF_onesec\2noth_ue.xls'];
fnm = [DATAPATH 'Entry3b_GG\Stat\DTF_onesec\stat_2noth_ue.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Noth\'];
mm = pwd;
cd(inpdir)
cd(inpplus)
files0 = dir(pwd);
files = {};
files_clus = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if isequal(ext,'.mat')
        files{end+1} = files0(i).name;
        files_clus{end+1} = [files0(i).name(1:end-7) 'CLUSTER'];
    end
end
sf = length(files);

% Statistics
NothNoburst = struct('name',{},'t_hypothesis_diff',{},'W_hypothesis_diff',{},...
    'KS_hypothesis_diff',{});
NothBurst = struct('name',{},'t_hypothesis_diff',{},'W_hypothesis_diff',{},...
    'KS_hypothesis_diff',{});
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
    PDC_diff_real = PDC_ue_real - PDC_eu_real;
    mPdiffr = mean(PDC_diff_real);
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
    PDC_diff_ctrl = PDC_ue_ctrl - PDC_eu_ctrl;
    mPdiffc = mean(PDC_diff_ctrl);
    mPuec = mean(PDC_ue_ctrl);
    alfa = 0.007;
    [KSh_diff,KSp_diff] = kstest2(PDC_diff_real,PDC_diff_ctrl,alfa,'smaller');
    [th_diff,tp_diff,ci_diff] = ttest(PDC_diff_real-PDC_diff_ctrl,0,alfa,'right');
    [Wp_diff,Wh_diff] = b_signrank2(PDC_diff_real,PDC_diff_ctrl,'alpha',alfa);
    KSh_diff = KSh_diff + 1 - 1;    % convert to numbers from logical values
    th_diff = th_diff + 1 - 1;
    Wh_diff = Wh_diff + 1 - 1;
    if mPdiffc > mPdiffr
        Wh_diff = 0;
        Wp_diff = -1;
    end
    if ClusterNumber == 0
        NothNoburst(end+1).name = files{o}(1:end-8);
        NothNoburst(end).t_hypothesis_diff = th_diff;
        NothNoburst(end).t_significance_diff = tp_diff;
        NothNoburst(end).W_hypothesis_diff = Wh_diff;
        NothNoburst(end).W_significance_diff = Wp_diff;
        NothNoburst(end).KS_hypothesis_diff = KSh_diff;
        NothNoburst(end).KS_significance_diff = KSp_diff;
        NothNoburst(end).mean_diff_real = mPdiffr;
        NothNoburst(end).mean_diff_ctrl = mPdiffc;
        NN{end+1,1} = th_diff;
        NN{end,2} = Wh_diff;
        NN{end,3} = KSh_diff;
        NN{end,4} = [];
        NN{end,5} = mPdiffr > mPdiffc;
        NN{end,5} = NN{end,5} + 1 - 1;
        NNname{end+1} = files{o}(1:end-8);
    else
        NothBurst(end+1).name = files{o}(1:end-8);
        NothBurst(end).t_hypothesis_diff = th_diff;
        NothBurst(end).t_significance_diff = tp_diff;
        NothBurst(end).W_hypothesis_diff = Wh_diff;
        NothBurst(end).W_significance_diff = Wp_diff;
        NothBurst(end).KS_hypothesis_diff = KSh_diff;
        NothBurst(end).KS_significance_diff = KSp_diff;
        NothBurst(end).mean_diff_real = mPdiffr;
        NothBurst(end).mean_diff_ctrl = mPdiffc;
        NB{end+1,1} = th_diff;
        NB{end,2} = Wh_diff;
        NB{end,3} = KSh_diff;
        NB{end,4} = [];
        NB{end,5} = mPdiffr > mPdiffc;
        NB{end,5} = NB{end,5} + 1 - 1;;
        NBname{end+1} = files{o}(1:end-8);
    end
end
alfa

% Save
str = [{'t_hypothesis_diff'} {'W_hypothesis_diff'} {'KS_hypothesis_diff'}...
    {' '} {'diff_r>c'}];
xlswrite(fn,str,'NothNoburst','B1');
xlswrite(fn,NNname','NothNoburst','A2');
xlswrite(fn,NN,'NothNoburst','B2');
xlswrite(fn,str,'NothBurst','B1');
xlswrite(fn,NBname','NothBurst','A2');
xlswrite(fn,NB,'NothBurst','B2');

save(fnm,'NothNoburst','NothBurst','NN','NB')
cd(mm)