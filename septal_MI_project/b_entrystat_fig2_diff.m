function [Stats, ST] = b_entrystat_fig2_diff
%ENTRYSTAT_FIG2_DIFF   Calculates statistics for entropy data.
%   ENTRYSTAT_FIG2_DIFF runs on the output of ENTRYRUN_THETA or ENTRYRUN3. 
%   Edit code to specify exact directories!
%
%   It runs on entropy results, and addresses the question whether HC-MS
%   and MS-HC information flow is larger under theta and under non-theta.
%   It calculates Wilcoxon signrank-test and saves results.
%
%   See also ENTRYSTAT5 and 2ENTRYSTAT_THETA4.

%Directories
global DATAPATH
inpdir_theta = [DATAPATH 'Entry3c\Theta\entropy\power\line\windowsize1000_overlap1\real\'];
inpdir_noth = [DATAPATH 'Entry3c\Noth\entropy\power\line\windowsize1000_overlap1\real\'];
ctrldir_theta = [DATAPATH 'Entry3c\Theta\entropy\power\line\windowsize1000_overlap1\control\'];
ctrldir_noth = [DATAPATH 'Entry3c\Noth\entropy\power\line\windowsize1000_overlap1\control\'];
fn = [DATAPATH 'Entrystat_fig2\newall_diff.xls'];
mm = pwd;

[files_theta, files_short_theta] = genfiles(inpdir_theta);
[files_noth, files_short_noth] = genfiles(inpdir_noth);
files_short = intersect(files_short_theta, files_short_noth);
sf = length(files_short);

% Statistics
Stats = struct([]);
ST = {};
STname = {};
for o = 1:sf
    cd(inpdir_theta)
    inx = find(strcmp(files_short{o},files_short_theta));
    ftheta = files_theta{inx};
    load(ftheta)
    ENT_eu_real_theta = aUxy;
    ENT_ue_real_theta = aUyx;
    mEeur_theta = median(aUxy);
    mEuer_theta = median(aUyx);
    mEdiff_theta = median(aUxy-aUyx);
    stEdiff_theta = std(aUxy-aUyx);

    cd(inpdir_noth)
    inx = find(strcmp(files_short{o},files_short_noth));
    fnoth = files_noth{inx};
    load(fnoth)
    ENT_eu_real_noth = aUxy;
    ENT_ue_real_noth = aUyx;
    mEeur_noth = median(aUxy);
    mEuer_noth = median(aUyx);
    mEdiff_noth = median(aUxy-aUyx);
    stEdiff_noth = std(aUxy-aUyx);
    
    cd(ctrldir_theta)
    load(ftheta)
    ENT_eu_ctrl_theta = aUxy;
    ENT_ue_ctrl_theta = aUyx;
    mEeur_theta_ctrl = median(aUxy);
    mEuer_theta_ctrl = median(aUyx);
    mEdiff_theta_ctrl = median(aUxy-aUyx);
    stEdiff_theta_ctrl = std(aUxy-aUyx);

    cd(ctrldir_noth)
    load(fnoth)
    ENT_eu_ctrl_noth = aUxy;
    ENT_ue_ctrl_noth = aUyx;
    mEeur_noth_ctrl = median(aUxy);
    mEuer_noth_ctrl = median(aUyx);
    mEdiff_noth_ctrl = median(aUxy-aUyx);
    stEdiff_noth_ctrl = std(aUxy-aUyx);
    
    if mEdiff_theta < 0
        Dd_theta = -mEdiff_theta + mEdiff_theta_ctrl;
    else
        Dd_theta = mEdiff_theta - mEdiff_theta_ctrl;
    end
    if mEdiff_noth < 0
        Dd_noth = -mEdiff_noth + mEdiff_noth_ctrl;
    else
        Dd_noth = mEdiff_noth - mEdiff_noth_ctrl;
    end

    if ~isequal(ftheta(1:6),fnoth(1:6))
        error('Technical error 60')
    end
    
    ST{end+1,1} = Dd_theta;
    ST{end,2} = Dd_noth;
    STname{end+1} = files_short{o};
end

% Save
str = [{'theta diff - theta diff ctrl'} {'noth diff - noth diff ctrl'}];
xlswrite(fn,str,'all','B1');
xlswrite(fn,STname','all','A2');
xlswrite(fn,ST,'all','B2');



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