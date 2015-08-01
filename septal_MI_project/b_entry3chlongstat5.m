function b_entry3chlongstat5
%ENTRY3CHLONGSTAT5   Calculates statistics for whole registration uncert. coeff.
%   ENTRY3CHLONGSTAT5 calculates and saves t-test, Wilcoxon signrank-test and
%   Kolmogorov-Smirnov test for whole registration uncertainty coefficients
%   (output of ENTRYRUN4). Edit code to modify directories!
%
%   Unlike ENTRYLONGSTAT, ENTRYLONGSTAT5 analyses segments with at least
%   100 spikes.
%
%   This program is for unit-unit analysis!
%
%   See also ENTRYLONGSTAT and ENTRYRUN4.

%   real_eu-real_ue>?ctrl_eu-ctrl_ue

% Directories
global DATAPATH
global DATADIR2
datadir = 'f:\raw_data\3ch_discriminated';
inpdir = [DATAPATH 'Entry_whole_3ch_Lu\entropy\power\line\real\temp\'];
inpdir2 = [DATAPATH 'Entry_whole_3ch_Lu\entropy\power\line\control\temp\'];
thetadir = [DATAPATH 'Wavelet_3ch\theta_segments\'];
nothdir = [DATAPATH 'Wavelet_3ch\nontheta_segments\'];
sharpdir = [DATAPATH 'Wavelet_3ch\sharpwave_segments\'];
resdir = [DATAPATH 'Entry_whole_3ch_Lu\Stat\'];
fnx = [resdir '\2longstat_eu.xls'];
fnm = [resdir '\2longstat_eu.mat'];
mm0 = pwd;

% Filelist
[files files_short] = b_filelist(inpdir);
files_short2 = unique(files_short);
sf = length(files_short2);
[datalist dlist] = flist2(datadir);

% Initialization
Theta = struct([]);
T = {};
Tname = {};
Tind1 = {};
Tind2 = {};
Noth = struct([]);
N = {};
Nname = {};
Nind1 = {};
Nind2 = {};
sr = 10000;     % sampling rate

% Progress indicator
wb = waitbar(0,'Running ENTRY3CHLONGSTAT5...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Load
for o = 1:sf
    fname = [files_short2{o} '_ENTROPYline.mat'];
    ff = fullfile(inpdir,fname);
    load(ff)
    aUxy_r = aUxy;
    aUyx_r = aUyx;
    ff = fullfile(inpdir2,fname);
    load(ff)
    aUxy_c = aUxy;
    aUyx_c = aUyx;

    [ftheta,ftheta_short] = flist(thetadir,1);
    [fnoth,fnoth_short] = flist(nothdir,2);
    [fsharp,fsharp_short] = flist(sharpdir,3);
    mm = pwd;
    cd(thetadir)
    inx = find(strcmp(fname(1:6),ftheta_short));
    load(ftheta(inx).name);
    cd(nothdir)
    inx = find(strcmp(fname(1:6),fnoth_short));
    load(fnoth(inx).name);
    cd(sharpdir)
    inx = find(strcmp(fname(1:6),fsharp_short));
    load(fsharp(inx).name);
    
    cd(datadir)
    idx = find(strcmp(files_short2{o},dlist));
    fn = datalist(idx(1)).name;
    load(fn)
    vdisc1 = vdisc;
    fn = datalist(idx(2)).name;
    load(fn)
    vdisc2 = vdisc;
    cd(mm)
    clear vdisc
    
% Statistics for theta segments
    if ~isempty(ThetaSegments)
        thetaseg = short_killer(min(ThetaSegments,length(aUxy)*sr),vdisc1);
        if ~isempty(thetaseg)
            thetaseg = short_killer(min(thetaseg,length(aUxy)*sr),vdisc2);
        end
    end
    for k1 = 1:size(thetaseg,2)
        ind1 = ceil(thetaseg(1,k1)/sr);
        ind2 = floor(thetaseg(2,k1)/sr);
        ENT_eu_real = aUxy_r(ind1:ind2);
        ENT_ue_real = aUyx_r(ind1:ind2);
        ENT_diff_real = ENT_eu_real - ENT_ue_real;
        mEeur = mean(ENT_eu_real);
        mEuer = mean(ENT_ue_real);
        mEdiffr = mean(ENT_diff_real);
        ENT_eu_ctrl = aUxy_c(ind1:ind2);
        ENT_ue_ctrl = aUyx_c(ind1:ind2);
        ENT_diff_ctrl = ENT_eu_ctrl - ENT_ue_ctrl;
        mEeuc = mean(ENT_eu_ctrl);
        mEuec = mean(ENT_ue_ctrl);
        mEdiffc = mean(ENT_diff_ctrl);
        alfa = 0.007;
        [KSh_diff,KSp_diff] = kstest2(ENT_diff_real,ENT_diff_ctrl,alfa,'smaller');
        [th_diff,tp_diff,ci_diff] = ttest(ENT_diff_real-ENT_diff_ctrl,0,alfa,'right');
        [Wp_diff,Wh_diff] = b_signrank2(ENT_diff_real,ENT_diff_ctrl,'alpha',alfa);
        KSh_diff = double(KSh_diff);
        th_diff = double(th_diff);
        Wh_diff = double(Wh_diff);
        if mEdiffc > mEdiffr
            Wh_diff = 0;
            Wp_diff = -1;
        end
        Theta(end+1).name = files_short2{o};
        Theta(end).ind1 = ind1;
        Theta(end).ind2 = ind2;
        Theta(end).t_hypothesis_diff = th_diff;
        Theta(end).t_significance_diff = tp_diff;
        Theta(end).W_hypothesis_diff = Wh_diff;
        Theta(end).W_significance_diff = Wp_diff;
        Theta(end).KS_hypothesis_diff = KSh_diff;
        Theta(end).KS_significance_diff = KSp_diff;
        Theta(end).mean_diff_real = mEdiffr;
        Theta(end).mean_diff_ctrl = mEdiffc;
        T{end+1,1} = th_diff;
        T{end,2} = Wh_diff;
        T{end,3} = KSh_diff;
        T{end,4} = [];
        T{end,5} = mEdiffr > mEdiffc;
        T{end,5} = double(T{end,5});
        Tname{end+1} = files_short2{o};
        Tind1{end+1} = ind1;
        Tind2{end+1} = ind2;
    end
    
% Statistics for non-theta segments
    if ~isempty(NonThetaSegments)
        nothseg = short_killer(min(NonThetaSegments,length(aUxy)*sr),vdisc1);
        if ~isempty(nothseg)
            nothseg = short_killer(min(nothseg,length(aUxy)*sr),vdisc2);
        end
    end
    for k2 = 1:size(nothseg,2)
        ind1 = ceil(nothseg(1,k2)/sr);
        ind2 = floor(nothseg(2,k2)/sr);
        ENT_eu_real = aUxy_r(ind1:ind2);
        ENT_ue_real = aUyx_r(ind1:ind2);
        ENT_diff_real = ENT_eu_real - ENT_ue_real;
        mEeur = mean(ENT_eu_real);
        mEuer = mean(ENT_ue_real);
        mEdiffr = mean(ENT_diff_real);
        ENT_eu_ctrl = aUxy_c(ind1:ind2);
        ENT_ue_ctrl = aUyx_c(ind1:ind2);
        ENT_diff_ctrl = ENT_eu_ctrl - ENT_ue_ctrl;
        mEeuc = mean(ENT_eu_ctrl);
        mEuec = mean(ENT_ue_ctrl);
        mEdiffc = mean(ENT_diff_ctrl);
        alfa = 0.007;
        [KSh_diff,KSp_diff] = kstest2(ENT_diff_real,ENT_diff_ctrl,alfa,'smaller');
        [th_diff,tp_diff,ci_diff] = ttest(ENT_diff_real-ENT_diff_ctrl,0,alfa,'right');
        [Wp_diff,Wh_diff] = b_signrank2(ENT_diff_real,ENT_diff_ctrl,'alpha',alfa);
        KSh_diff = double(KSh_diff);
        th_diff = double(th_diff);
        Wh_diff = double(Wh_diff);
        if mEdiffc > mEdiffr
            Wh_diff = 0;
            Wp_diff = -1;
        end
        Noth(end+1).name = files_short2{o};
        Noth(end).ind1 = ind1;
        Noth(end).ind2 = ind2;
        Noth(end).t_hypothesis_diff = th_diff;
        Noth(end).t_significance_diff = tp_diff;
        Noth(end).W_hypothesis_diff = Wh_diff;
        Noth(end).W_significance_diff = Wp_diff;
        Noth(end).KS_hypothesis_diff = KSh_diff;
        Noth(end).KS_significance_diff = KSp_diff;
        Noth(end).mean_diff_real = mEdiffr;
        Noth(end).mean_diff_ctrl = mEdiffc;
        N{end+1,1} = th_diff;
        N{end,2} = Wh_diff;
        N{end,3} = KSh_diff;
        N{end,4} = [];
        N{end,5} = mEdiffr > mEdiffc;
        N{end,5} = double(N{end,5});
        Nname{end+1} = files_short2{o};
        Nind1{end+1} = ind1;
        Nind2{end+1} = ind2;
    end
    waitbar(o/sf)
end
alfa
close(wb)

% Save
str = [{' '} {' '} {'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'diff_r>?c'}];
xlswrite(fnx,str,'Theta','B1');
xlswrite(fnx,Tname','Theta','A2');
xlswrite(fnx,Tind1','Theta','B2');
xlswrite(fnx,Tind2','Theta','C2');
xlswrite(fnx,T,'Theta','D2');

xlswrite(fnx,str,'Noth','B1');
xlswrite(fnx,Nname','Noth','A2');
xlswrite(fnx,Nind1','Noth','B2');
xlswrite(fnx,Nind2','Noth','C2');
xlswrite(fnx,N,'Noth','D2');

save(fnm,'Theta','T')
save(fnm,'Noth','N')
cd(mm0)



% -------------------------------------------------------------------------
function [files2,files2_short] = flist(inpdir,k);

switch k
    case 1
        c1 = 16;
        c2 = 21;
    case 2
        c1 = 19;
        c2 = 24;
    case 3
        c1 = 20;
        c2 = 25;
end

files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name(c1:min(c2,length(files(i).name)));
    end
end
files2 = files2(2:end);



% -------------------------------------------------------------------------
function [files2, files2_short] = flist2(inpdir)

files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:3) '_n' files(i).name(6)];
    end
end
files2 = files2(2:end);



% ----------------------------------------------------------------------------------
function segments = short_killer(segments,vdisc)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<50000);         % leaving segments shorter than 5 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];

% Skip segments with less than 100 spikes
int = segments;
int1 = int(1,:);
int2 = int(2,:);
fd = [];
for k = 1:length(int1)
    lvd = length(find(vdisc>int1(k)&vdisc<int2(k)));
    if lvd < 100
        fd(end+1) = k;
    end
end
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];