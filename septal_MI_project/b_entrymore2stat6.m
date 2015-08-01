function [Stats_eu,STeu] = b_entrymore2stat6
%ENTRYMORE2STAT6   Calculates statistics for entropy data.
%   ENTRYMORE2STAT6 runs on the output of ENTRYRUN_THETA or ENTRYRUN3C. Edit
%   code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS
%   and MS-HC information flow is larger under theta or not. It calculates
%   and saves 3 statistics: t-test, Wilcoxon ranksum-test (Mann-Whitney
%   U-test) and Kolmogorov-Smirnov test.
%
%   It analysis mutual information values instead of uncertainty coefficients.
%
%   See also ENTRYSTAT5.

% entropy
% theta and noth
% growing?

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'Entry3b_M\Theta\entropy\power\line\windowsize1000_overlap1\real\'];
inpdir_noth = [DATAPATH 'Entry3b_M\Noth\entropy\power\line\windowsize1000_overlap1\real\'];
fn = [DATAPATH 'Entry3b_M\Stat\power\2all_mutinf.xls'];
fnm = [DATAPATH 'Entry3b_M\Stat\power\stat_2all_mutinf.mat'];
mm = pwd;

[files_theta, files_short_theta] = genfiles(inpdir_theta);
[files_noth, files_short_noth] = genfiles(inpdir_noth);
files_short = intersect(files_short_theta, files_short_noth);
sf = length(files_short);

% Statistics
Stats_eu = struct([]);
STeu = {};
STname = {};
for o = 1:sf
    cd(inpdir_theta)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    ENT_eu_real_theta = aIxy;
    mEeur_theta = mean(aIxy);
    
    cd(inpdir_noth)
    inx = find(strcmp(files_short{o},files_short_noth));
    fnoth = files_noth{inx};
    load(fnoth)
    ENT_eu_real_noth = aIxy;
    mEeur_noth = mean(aIxy);
    
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
    
    STname{end+1} = files_short{o};
end
alfa

% Save
str = [{'t Ieu'} {'W Ieu'} {'KS Ieu'}...
    {' '} {'Ieu_theta>?Ieu_noth'}];
xlswrite(fn,str,'growing','B1');
xlswrite(fn,STname','growing','A2');
xlswrite(fn,STeu,'growing','B2');

save(fnm,'Stats_eu','STeu')
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