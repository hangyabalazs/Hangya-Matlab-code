function [Stats_eu, Stats_ue, STeu, STue] = b_entrystat6
%ENTRYSTAT6   Calculates statistics for entropy data.
%   ENTRYSTAT6 runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS
%   and MS-HC information flow is larger under theta or not. It calculates
%   and saves 3 statistics: t-test, Wilcoxon ranksum-test (Mann-Whitney
%   U-test) and Kolmogorov-Smirnov test.
%
%   See also ENTRYSTAT5.

% entropy
% theta and noth
% growing?

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'Entry3c_cont3_unitunit\Theta\entropy\phase\line\windowsize1000_overlap1\real\'];
inpdir_noth = [DATAPATH 'Entry3c_cont3_unitunit\Noth\entropy\phase\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3c_cont3_unitunit\Stat\phase\2all.xls'];
fnm = [DATAPATH 'Entry3c_cont3_unitunit\Stat\phase\stat_2all.mat'];
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
    ENT_eu_real_theta = aUxy;
    ENT_ue_real_theta = aUyx;
    mEeur_theta = mean(aUxy);
    mEuer_theta = mean(aUyx);

    cd(inpdir_noth)
    inx = find(strcmp(files_short{o},files_short_noth));
    fnoth = files_noth{inx};
    load(fnoth)
    ENT_eu_real_noth = aUxy;
    ENT_ue_real_noth = aUyx;
    mEeur_noth = mean(aUxy);
    mEuer_noth = mean(aUyx);

    if ~isequal(ftheta(1:6),fnoth(1:6))
        error('Technical error 40')
    end
    alfa = 0.007;
    if mEeur_theta > mEeur_noth
        [KSh,KSp] = kstest2(ENT_eu_real_theta,ENT_eu_real_noth,alfa,'smaller');
        [th,tp,ci] = ttest2(ENT_eu_real_theta,ENT_eu_real_noth,alfa,'right','unequal');
        [Wp,Wh] = b_ranksum2(ENT_eu_real_theta,ENT_eu_real_noth,'alpha',alfa);
    else
        [KSh,KSp] = kstest2(ENT_eu_real_theta,ENT_eu_real_noth,alfa,'larger');
        [th,tp,ci] = ttest2(ENT_eu_real_theta,ENT_eu_real_noth,alfa,'left','unequal');
        [Wp,Wh] = b_ranksum2(ENT_eu_real_theta,ENT_eu_real_noth,'alpha',alfa);
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
    Stats_eu(end).mean_eu_theta = mEeur_theta;
    Stats_eu(end).mean_eu_noth = mEeur_noth;
    STeu{end+1,1} = th;
    STeu{end,2} = Wh;
    STeu{end,3} = KSh;
    STeu{end,4} = [];
    STeu{end,5} = mEeur_theta > mEeur_noth;
    STeu{end,5} = STeu{end,5} + 1 - 1;
    STeu{end,6} = [];
    
    if mEuer_theta > mEuer_noth
        [KSh,KSp] = kstest2(ENT_ue_real_theta,ENT_ue_real_noth,alfa,'smaller');
        [th,tp,ci] = ttest2(ENT_ue_real_theta,ENT_ue_real_noth,alfa,'right','unequal');
        [Wp,Wh] = b_ranksum2(ENT_ue_real_theta,ENT_ue_real_noth,'alpha',alfa);
    else
        [KSh,KSp] = kstest2(ENT_ue_real_theta,ENT_ue_real_noth,alfa,'larger');
        [th,tp,ci] = ttest2(ENT_ue_real_theta,ENT_ue_real_noth,alfa,'left','unequal');
        [Wp,Wh] = b_ranksum2(ENT_ue_real_theta,ENT_ue_real_noth,'alpha',alfa);
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
    Stats_ue(end).mean_ue_theta = mEuer_theta;
    Stats_ue(end).mean_ue_noth = mEuer_noth;
    STue{end+1,1} = th;
    STue{end,2} = Wh;
    STue{end,3} = KSh;
    STue{end,4} = [];
    STue{end,5} = mEuer_theta > mEuer_noth;
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