function b_entryzshiftstatUE(k)
%ENTRYZSHIFTSTATUE   Calculates statistics for zshifted entropy data.
%   ENTRYZSHIFTSTATUE runs on the output of ENTRYZSHIFT_COMPARE. Edit
%   code to specify exact directories!
%
%   It runs on entropy and zshifted entropy results, and addresses the
%   question whether zshifted entropy is higher than entropy without zshift 
%   or not. It calculates and saves 3 statistics: t-test, Wilcoxon ranksum-
%   test (Mann-Whitney U-test) and Kolmogorov-Smirnov test.
%
%   ENTRYZSHIFTSTATUE(K) calculates with the first K values of entropy
%   data.
%
%   It calculates statistics for MS-HC direction!
%
%   See also ENTRYZSHIFT_COMPARE.

% Input argument check
error(nargchk(0,1,nargin))
if nargin < 1
    k = Inf;
end

%Directories
name = '';
global DATAPATH
inpdir = [DATAPATH 'Entry_zshift_onset\Compare\mat\'];
if isequal(k,Inf)
    fn = [DATAPATH 'Entry_zshift_onset\Stat\UEstat.xls'];
else
    fn = [DATAPATH 'Entry_zshift_onset\Stat\UEstat_' num2str(k) '.xls'];
end
mm = pwd;
cd(inpdir)
files = b_filelist(pwd);
sf = length(files);

% Statistics
T = {};
Tname = {};
for o = 1:sf
    cd(inpdir)
    load(files(o).name)
    lUE1 = length(UnitEeg1);
    lUE2 = length(UnitEeg2);
    mUE1 = mean(UnitEeg1(1:min(k,lUE1)));
    mUE2 = mean(UnitEeg2(1:min(k,lUE2)));
    if mUE2 > mUE1
        [KSh,KSp] = kstest2(UnitEeg2(1:min(k,lUE2)),UnitEeg1(1:min(k,lUE1)),0.05,'smaller');
        [th,tp,ci] = ttest2(UnitEeg2(1:min(k,lUE2)),UnitEeg1(1:min(k,lUE1)),[],'right','unequal');
        [Wp,Wh] = b_ranksum2(UnitEeg2(1:min(k,lUE2)),UnitEeg1(1:min(k,lUE1)));
    else
        [KSh,KSp] = kstest2(UnitEeg2(1:min(k,lUE2)),UnitEeg1(1:min(k,lUE1)),0.05,'larger');
        [th,tp,ci] = ttest2(UnitEeg2(1:min(k,lUE2)),UnitEeg1(1:min(k,lUE1)),[],'left','unequal');
        [Wp,Wh] = b_ranksum2(UnitEeg2(1:min(k,lUE2)),UnitEeg1(1:min(k,lUE1)));
    end
    KSh = KSh + 1 - 1;    % convert to numbers from logical values
    th = th + 1 - 1;
    Wh = Wh + 1 - 1;
    T{end+1,1} = th;
    T{end,2} = Wh;
    T{end,3} = KSh;
    T{end,4} = [];
    T{end,5} = mUE2 > mUE1;
    T{end,5} = T{end,5} + 1 - 1;
    Tname{end+1} = files(o).name(1:end-16);
end

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'UE2>?UE1'}];
xlswrite(fn,str,'EntryZshift','B1');
xlswrite(fn,Tname','EntryZshift','A2');
xlswrite(fn,T,'EntryZshift','B2');