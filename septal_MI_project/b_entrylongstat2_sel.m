function b_entrylongstat2_sel
%ENTRYLONGSTAT2_SEL   Calculates statistics for whole registration uncert. coeff.
%   ENTRYLONGSTAT2_SEL is a modified version of ENTRYLONGSTAT2 for long
%   segments that are "selected", i.e. a 175 sec long segment is selected
%   from the registration (not starting at the beginning of the recording).
%   See ENTRYLONGSTAT2 for details.
%
%   See also ENTRYLONGSTAT2.

% Directories
global DATAPATH
global DATADIR2
datadir = DATADIR2;
inpdir = [DATAPATH 'Entry_whole_selected\entropy\power\line\real\'];
thetadir = [DATAPATH 'Wavelet2\theta_segments\'];
nothdir = [DATAPATH 'Wavelet2\nontheta_segments\'];
sharpdir = [DATAPATH 'Wavelet2\sharpwave_segments\'];
resdir = [DATAPATH 'Entry_whole\Stat\'];
fnx = [resdir '\2longstat_sel.xls'];
fnm = [resdir '\2longstat_sel.mat'];
mm0 = pwd;

% Filelist
[files files_short] = b_filelist(inpdir);
files_short2 = unique(files_short);
sf = length(files_short2);
[datalist dlist] = b_filelist(datadir);

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
wb = waitbar(0,'Running ENTRYLONGSTAT...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Correction for selected segment
fst = 35;
lst = fst + 175;

% Load
for o = 4:4
    fname = [files_short2{o} '_ENTROPYline.mat'];
    ff = fullfile(inpdir,fname);
    load(ff)

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
    fn = datalist(idx).name;
    load(fn)
    cd(mm)
    
% Statistics for theta segments
    thetaseg = short_killer(min(max(ThetaSegments,fst*sr),lst*sr),vdisc);
    for k1 = 1:size(thetaseg,2)
        ind1 = - fst + ceil(thetaseg(1,k1)/sr);
        ind2 = - fst + floor(thetaseg(2,k1)/sr);
        ENT_eu_real = aUxy(ind1:ind2);
        ENT_ue_real = aUyx(ind1:ind2);
        mEeur = mean(ENT_eu_real);
        mEuer = mean(ENT_ue_real);
        if mEeur > mEuer
            [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,0.05,'smaller');
            [th,tp,ci] = ttest(ENT_eu_real-ENT_ue_real,0,[],'right');
            [Wp,Wh] = b_signrank2(ENT_eu_real-ENT_ue_real);
        else
            [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,0.05,'larger');
            [th,tp,ci] = ttest(ENT_eu_real-ENT_ue_real,0,[],'left');
            [Wp,Wh] = b_signrank2(ENT_ue_real-ENT_eu_real);
        end
        KSh = double(KSh);
        th = double(th);
        Wh = double(Wh);
        Theta(end+1).name = files_short2{o};
        Theta(end).ind1 = ind1 + fst;
        Theta(end).ind2 = ind2 + fst;
        Theta(end).t_hypothesis = th;
        Theta(end).t_significance = tp;
        Theta(end).W_hypothesis = Wh;
        Theta(end).W_significance = Wp;
        Theta(end).KS_hypothesis = KSh;
        Theta(end).KS_significance = KSp;
        Theta(end).mean_eu_real = mEeur;
        Theta(end).mean_ue_real = mEuer;
        T{end+1,1} = th;
        T{end,2} = Wh;
        T{end,3} = KSh;
        T{end,4} = [];
        T{end,5} = mEeur > mEuer;
        T{end,5} = double(T{end,5});
        Tname{end+1} = files_short2{o};
        Tind1{end+1} = ind1 + fst;
        Tind2{end+1} = ind2 + fst;
    end
    
% Statistics for non-theta segments
    nothseg = short_killer(min(max(NonThetaSegments,fst*sr),lst*sr),vdisc);
    for k2 = 1:size(nothseg,2)
        
        ind1 = - fst + ceil(nothseg(1,k2)/sr);
        ind2 = - fst + floor(nothseg(2,k2)/sr);
        ENT_eu_real = aUxy(ind1:ind2);
        ENT_ue_real = aUyx(ind1:ind2);
        mEeur = mean(ENT_eu_real);
        mEuer = mean(ENT_ue_real);
        if mEeur > mEuer
            [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,0.05,'smaller');
            [th,tp,ci] = ttest(ENT_eu_real-ENT_ue_real,0,[],'right');
            [Wp,Wh] = b_signrank2(ENT_eu_real-ENT_ue_real);
        else
            [KSh,KSp] = kstest2(ENT_eu_real,ENT_ue_real,0.05,'larger');
            [th,tp,ci] = ttest(ENT_eu_real-ENT_ue_real,0,[],'left');
            [Wp,Wh] = b_signrank2(ENT_ue_real-ENT_eu_real);
        end
        KSh = double(KSh);
        th = double(th);
        Wh = double(Wh);
        Noth(end+1).name = files_short2{o};
        Noth(end).ind1 = ind1 + fst;
        Noth(end).ind2 = ind2 + fst;
        Noth(end).t_hypothesis = th;
        Noth(end).t_significance = tp;
        Noth(end).W_hypothesis = Wh;
        Noth(end).W_significance = Wp;
        Noth(end).KS_hypothesis = KSh;
        Noth(end).KS_significance = KSp;
        Noth(end).mean_eu_real = mEeur;
        Noth(end).mean_ue_real = mEuer;
        N{end+1,1} = th;
        N{end,2} = Wh;
        N{end,3} = KSh;
        N{end,4} = [];
        N{end,5} = mEeur > mEuer;
        N{end,5} = double(N{end,5});
        Nname{end+1} = files_short2{o};
        Nind1{end+1} = ind1 + fst;
        Nind2{end+1} = ind2 + fst;
    end
    waitbar(o/sf)
end
close(wb)

% Save
str = [{'t_hypothesis'} {'W_hypothesis'} {'KS_hypothesis'}...
    {' '} {'eu>?ue'}];
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