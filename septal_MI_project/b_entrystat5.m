function [Stats_eu, Stats_ue, STeu, STue] = b_entrystat5
%ENTRYSTAT5   Calculates statistics for entropy data.
%   ENTRYSTAT5 runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact directories!
%
%   It runs on DTF results, and addresses the question whether HC-MS
%   and MS-HC information flow is larger under theta or not. It calculates
%   and saves 3 statistics: t-test, Wilcoxon ranksum-test (Mann-Whitney
%   U-test) and Kolmogorov-Smirnov test.
%
%   See also ENTRYSTAT6.

% DTF
% theta and noth
% growing?

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'Entry3_cont2\Theta\DTF\onesec\real\'];
inpdir_noth = [DATAPATH 'Entry3_cont2\Noth\DTF\onesec\real\'];
fn = [DATAPATH 'Entry3_cont2\Stat\DTF_onesec\2all.xls'];
fnm = [DATAPATH 'Entry3_cont2\Stat\DTF_onesec\stat_2all.mat'];
mm = pwd;

[files_theta, files_short_theta] = genfiles(inpdir_theta);
[files_noth, files_short_noth] = genfiles(inpdir_noth);
files_short = intersect(files_short_theta, files_short_noth);
sf = length(files_short);

% Statistics
Stats_eu = struct([]);
Stats_ue = struct([]);
STeu = {};
STue = {};
STname = {};
for o = 1:sf
    cd(inpdir_theta)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    PDC_eu_real_theta = mean(PDC_eu(3:6,:));
    PDC_ue_real_theta = mean(PDC_ue(3:6,:));
    mPeur_theta = mean(PDC_eu_real_theta);
    mPuer_theta = mean(PDC_ue_real_theta);

    cd(inpdir_noth)
    inx = find(strcmp(files_short{o},files_short_noth));
    fnoth = files_noth{inx};
    load(fnoth)
    PDC_eu_real_noth = mean(PDC_eu(3:6,:));
    PDC_ue_real_noth = mean(PDC_ue(3:6,:));
    mPeur_noth = mean(PDC_eu_real_noth);
    mPuer_noth = mean(PDC_ue_real_noth);

    if ~isequal(ftheta(1:6),fnoth(1:6))
        error('Technical error 40')
    end
    alfa = 0.007;
    if mPeur_theta > mPeur_noth
        [KSh,KSp] = kstest2(PDC_eu_real_theta,PDC_eu_real_noth,alfa,'smaller');
        [th,tp,ci] = ttest2(PDC_eu_real_theta,PDC_eu_real_noth,alfa,'right','unequal');
        [Wp,Wh] = b_ranksum2(PDC_eu_real_theta,PDC_eu_real_noth,'alpha',alfa);
    else
        [KSh,KSp] = kstest2(PDC_eu_real_theta,PDC_eu_real_noth,alfa,'larger');
        [th,tp,ci] = ttest2(PDC_eu_real_theta,PDC_eu_real_noth,alfa,'left','unequal');
        [Wp,Wh] = b_ranksum2(PDC_eu_real_theta,PDC_eu_real_noth,'alpha',alfa);
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    Stats_eu(end+1).name = files_short{o};
    Stats_eu(end).t_hypothesis = th;
    Stats_eu(end).t_significance = tp;
    Stats_eu(end).W_hypothesis = Wh;
    Stats_eu(end).W_significance = Wp;
    Stats_eu(end).KS_hypothesis = KSh;
    Stats_eu(end).KS_significance = KSp;
    Stats_eu(end).mean_eu_theta = mPeur_theta;
    Stats_eu(end).mean_eu_noth = mPeur_noth;
    STeu{end+1,1} = th;
    STeu{end,2} = Wh;
    STeu{end,3} = KSh;
    STeu{end,4} = [];
    STeu{end,5} = mPeur_theta > mPeur_noth;
    STeu{end,5} = STeu{end,5} + 1 - 1;
    STeu{end,6} = [];

    if mPuer_theta > mPuer_noth
        [KSh,KSp] = kstest2(PDC_ue_real_theta,PDC_ue_real_noth,alfa,'smaller');
        [th,tp,ci] = ttest2(PDC_ue_real_theta,PDC_ue_real_noth,alfa,'right','unequal');
        [Wp,Wh] = b_ranksum2(PDC_ue_real_theta,PDC_ue_real_noth,'alpha',alfa);
    else
        [KSh,KSp] = kstest2(PDC_ue_real_theta,PDC_ue_real_noth,alfa,'larger');
        [th,tp,ci] = ttest2(PDC_ue_real_theta,PDC_ue_real_noth,alfa,'left','unequal');
        [Wp,Wh] = b_ranksum2(PDC_ue_real_theta,PDC_ue_real_noth,'alpha',alfa);
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    Stats_ue(end+1).name = files_short{o};
    Stats_ue(end).t_hypothesis = th;
    Stats_ue(end).t_significance = tp;
    Stats_ue(end).W_hypothesis = Wh;
    Stats_ue(end).W_significance = Wp;
    Stats_ue(end).KS_hypothesis = KSh;
    Stats_ue(end).KS_significance = KSp;
    Stats_ue(end).mean_ue_theta = mPuer_theta;
    Stats_ue(end).mean_ue_noth = mPuer_noth;
    STue{end+1,1} = th;
    STue{end,2} = Wh;
    STue{end,3} = KSh;
    STue{end,4} = [];
    STue{end,5} = mPuer_theta > mPuer_noth;
    STue{end,5} = STue{end,5} + 1 - 1;
    STname{end+1} = files_short{o};
end
alfa

% Save
str = [{'t eu'} {'W eu'} {'KS eu'}...
    {' '} {'eu_theta>?eu_noth'} {' '}...
    {'t ue'} {'W ue'} {'KS ue'}...
    {' '} {'ue_theta>?ue_noth'}];
xlswrite(fn,str,'growing','B1');
xlswrite(fn,STname','growing','A2');
xlswrite(fn,[STeu STue],'growing','B2');

save(fnm,'Stats_eu','Stats_ue','STeu','STue')
cd(mm)



% -------------------------------------------------------------------------
function [files, files_short] = genfiles(inpdir)

mm = pwd;
cd(inpdir)
files0 = dir(pwd);
files = {};
files_short = {};
for i = 1:length(files0)
    [pth name ext] = fileparts(files0(i).name);
    if isequal(ext,'.mat')
        files{end+1} = files0(i).name;
        files_short{end+1} = files0(i).name(1:6);
    end
end
cd(mm)