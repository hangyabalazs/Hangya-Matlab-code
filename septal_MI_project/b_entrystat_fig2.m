function [Stats, ST] = b_entrystat_fig2
%ENTRYSTAT_FIG2   Calculates statistics for entropy data.
%   ENTRYSTAT_FIG2 runs on the output of ENTRYRUN_THETA or ENTRYRUN3. Edit
%   code to specify exact directories!
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
fn = [DATAPATH 'Entrystat_fig2\newall.xls'];
fnm = [DATAPATH 'Entrystat_fig2\stat_newall.mat'];
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

    if ~isequal(ftheta(1:6),fnoth(1:6))
        error('Technical error 60')
    end
    
    if mEeur_theta > mEuer_theta
        [Wp_theta,Wh_theta] = b_signrank2(ENT_eu_real_theta-ENT_ue_real_theta);
    else
        [Wp_theta,Wh_theta] = b_signrank2(ENT_ue_real_theta-ENT_eu_real_theta);
    end
    if mEeur_noth > mEuer_noth
        [Wp_noth,Wh_noth] = b_signrank2(ENT_eu_real_noth-ENT_ue_real_noth);
    else
        [Wp_noth,Wh_noth] = b_signrank2(ENT_ue_real_noth-ENT_eu_real_noth);
    end
    Wh_theta = double(Wh_theta);    % convert to numbers from logical values
    Wh_noth = double(Wh_noth);
    Stats(end+1).name = files_short{o};
    Stats(end).mean_diff_theta = mEdiff_theta;
    Stats(end).st_diff_theta = stEdiff_theta;
    Stats(end).mean_diff_noth = mEdiff_noth;
    Stats(end).st_diff_noth = stEdiff_noth;
    Stats(end).mean_eu_noth = mEeur_noth;
    Stats(end).Wh_theta = Wh_theta;
    Stats(end).Wp_theta = Wp_theta;
    Stats(end).Wh_noth = Wh_noth;
    Stats(end).Wp_noth = Wp_noth;
    ST{end+1,1} = mEdiff_theta;
    ST{end,2} = stEdiff_theta;
    ST{end,3} = mEdiff_noth;
    ST{end,4} = stEdiff_noth;
    ST{end,5} = mEeur_theta > mEuer_theta;
    ST{end,5} = double(ST{end,5});
    ST{end,6} = mEeur_noth > mEuer_noth;
    ST{end,6} = double(ST{end,6});
    ST{end,7} = Wh_theta;
    ST{end,8} = Wp_theta;
    ST{end,9} = Wh_noth;
    ST{end,10} = Wp_noth;
    
    STname{end+1} = files_short{o};
end

% Save
str = [{'theta mean'} {'theta std'} {'noth mean'} {'noth std'}...
    {'eu_theta>?ue_theta'} {'eu_noth>?ue_noth'}...
    {'Wh theta'} {'Wp theta'} {'Wh noth'} {'Wp noth'}];
xlswrite(fn,str,'all','B1');
xlswrite(fn,STname','all','A2');
xlswrite(fn,ST,'all','B2');

save(fnm,'Stats','ST')



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