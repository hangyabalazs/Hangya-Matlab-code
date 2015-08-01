function [ThetaNoburst,ThetaBurst,TN,TB] = b_entrystat3_theta
%ENTRYSTAT3_THETA   Calculates statistics for entropy data (theta segments).
%   ENTRYSTAT3_THETA runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact direcctories!
%
%   It runs on PDC results, and addresses the question whether HC-MS 
%   information flow is larger than MS-HC or not. It calculates and saves
%   3 statistics: t-test, Wilcoxon ranksum-test (Mann-Whitney U-test) and
%   Kolmogorov-Smirnov test.
%
%   See also ENTRYSTAT3_NOTH and ENTRYSTAT4_THETA.

% DTF
% theta
% eu>?ue

%Directories
global DATAPATH
inpdir = [DATAPATH 'Entry2_unitunit\'];
inpplus = ['Theta\dtf\halfsec\real\'];
fn = [DATAPATH 'Entry2_unitunit\Stat\DTF_halfsec\2theta.xls'];
fnm = [DATAPATH 'Entry2_unitunit\Stat\DTF_halfsec\stat_2theta.mat'];
clusdir = [DATAPATH 'Burst\Cluster\Theta_3ch\'];
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
ThetaNoburst = struct('name',{},'t_hypothesis',{},'W_hypothesis',{},...
    'KS_hypothesis',{});
ThetaBurst = struct('name',{},'t_hypothesis',{},'W_hypothesis',{},...
    'KS_hypothesis',{});
TN = {};
TB = {};
TNname = {};
TBname = {};
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
    if mPeur > mPuer
        [KSh,KSp] = kstest2(PDC_eu_real,PDC_ue_real,0.05,'smaller');
        [th,tp,ci] = ttest2(PDC_eu_real,PDC_ue_real,[],'right','unequal');
        [Wp,Wh] = b_ranksum2(PDC_eu_real,PDC_ue_real);
    else
        [KSh,KSp] = kstest2(PDC_eu_real,PDC_ue_real,0.05,'larger');
        [th,tp,ci] = ttest2(PDC_eu_real,PDC_ue_real,[],'left','unequal');
        [Wp,Wh] = b_ranksum2(PDC_ue_real,PDC_eu_real);
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    if ClusterNumber == 0
        ThetaNoburst(end+1).name = files{o}(1:end-8);
        ThetaNoburst(end).t_hypothesis = th;
        ThetaNoburst(end).t_significance = tp;
        ThetaNoburst(end).W_hypothesis = Wh;
        ThetaNoburst(end).W_significance = Wp;
        ThetaNoburst(end).KS_hypothesis = KSh;
        ThetaNoburst(end).KS_significance = KSp;
        ThetaNoburst(end).mean_eu_real = mPeur;
        ThetaNoburst(end).mean_ue_real = mPuer;
        TN{end+1,1} = th;
        TN{end,2} = Wh;
        TN{end,3} = KSh;
        TN{end,4} = [];
        TN{end,5} = mPeur > mPuer;
        TN{end,5} = TN{end,5} + 1 - 1;
        TNname{end+1} = files{o}(1:end-8);
    else
        ThetaBurst(end+1).name = files{o}(1:end-8);
        ThetaBurst(end).t_hypothesis = th;
        ThetaBurst(end).t_significance = tp;
        ThetaBurst(end).W_hypothesis = Wh;
        ThetaBurst(end).W_significance = Wp;
        ThetaBurst(end).KS_hypothesis = KSh;
        ThetaBurst(end).KS_significance = KSp;
        ThetaBurst(end).mean_eu_real = mPeur;
        ThetaBurst(end).mean_ue_real = mPuer;
        TB{end+1,1} = th;
        TB{end,2} = Wh;
        TB{end,3} = KSh;
        TB{end,4} = [];
        TB{end,5} = mPeur > mPuer;
        TB{end,5} = TB{end,5} + 1 - 1;
        TBname{end+1} = files{o}(1:end-8);
    end
end

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'eu>?ue'}];
xlswrite(fn,str,'ThetaNoburst','B1');
xlswrite(fn,TNname','ThetaNoburst','A2');
xlswrite(fn,TN,'ThetaNoburst','B2');
xlswrite(fn,str,'ThetaBurst','B1');
xlswrite(fn,TBname','ThetaBurst','A2');
xlswrite(fn,TB,'ThetaBurst','B2');

save(fnm,'ThetaNoburst','ThetaBurst','TN','TB')