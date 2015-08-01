function arestrictstat_burstphase3
%ARESTRICTSTAT_BURSTPHASE3   Summary phase statistics of ordinal intraburst spikes in bounded frequency band.
%   ARESTRICTSTAT_BURSTPHASE3 calculates summary phase statistics of
%   ordinal intraburst spikes for Po, VPM, PoVPM, LD and nRT data. It plots
%   and saves result figures displaying differences of mean ordinal spike
%   phases. A structure ('dfstat') containing the number and proportion of
%   positive differences as well as the overall mean difference is also
%   saved. Boxplot and bargraphs of pooled differences are saved and their
%   relation to zero median is tested using Mann-Whitney U-test (note, that
%   the issue of concidering a narrow range circular sample as linear
%   should be discussed here). Test results appear on the boxplot as well
%   as in the 'dfstat' structure. Edit code to modify
%   input and output directories!
%
%   See also ARESTRICTSTAT_BURST and ARESTRICTSTAT_PHASE3.

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_stand\'];   % phase analysis data
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_stand\'];
mm = pwd;
dr = dir(inpdir);

% Main
MRL = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
MN = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
DF = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
for o = 3:length(dr)
    inpadd = dr(o).name;    % load burst data
    cd(inpdir)
    cd(inpadd)
    cd('bas')
    ddr = dir(pwd);
    fn = ddr(end).name(1:end-4);
    cmps = strread(fn,'%s','delimiter','_');
    fname = [cmps{1} '_' cmps{2}];
    tstr = [cmps{1} ' ' cmps{2}];
    ff1 = [inpdir2 fn '_PHASE.mat'];
    ff2 = [inpdir2 fn '_BURSTPHASE.mat'];
    try
        load(ff1)
        load(ff2)
    catch
        lasterr
        continue
    end
    ff = [tabledir 'tablazat_Balazsnak'];   % load position data
    [tbl0 tbl] = xlsread(ff);
    inx = find(strcmp({tbl{:,1}},fname));
    loc = tbl{inx,3};
    
    ibspno = H1ibspno;
    mibs = max(ibspno);
    edges = [-180:20:180];
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    ftm = zeros(1,mibs-1);
    mn_rad = zeros(1,mibs-1);
    mn = zeros(1,mibs-1);
    dfmn = zeros(1,mibs-2);
    mrl = zeros(1,mibs-1);
    dsc = 10;
    for k = 2:mibs
        str = ['let = find(ibspno==' num2str(k) ');'];
        eval(str)
        lele = length(let);
        ng = aang_fs(let) / 180 * pi;
        ftm(k-1) = sum(exp(1).^(i*ng)) / lele;
        if lele > dsc
            mn_rad(k-1) = angle(ftm(k-1));
            mn(k-1) = angle(ftm(k-1)) * 180 / pi;
            mrl(k-1) = abs(ftm(k-1));
            nt = eval(['length(MRL.' loc ');']);
            if k - 1 > nt
                ntt = 0;
            else
                ntt = eval(['length(MRL.' loc '{k-1});']);
            end
            eval(['MRL.' loc '{k-1}(ntt+1) = mrl(k-1);']);
            eval(['MN.' loc '{k-1}(ntt+1) = mn(k-1);']);
        else
            mn_rad(k-1) = NaN;
            mn(k-1) = NaN;
            mrl(k-1) = NaN;
        end
        if k > 2
            dfmn(k-2) = circdiff(mn(k-1),mn(k-2));
        end
    end
    for k = 2:mibs-1
        nt = eval(['length(DF.' loc ');']);
        if k - 1 > nt
            ntt = 0;
        else
            ntt = eval(['length(DF.' loc '{k-1});']);
        end
        if ~isnan(dfmn(k-1))
            eval(['DF.' loc '{k-1}(ntt+1) = dfmn(k-1);']);
        end
    end
end

% Saving statistics
DFPo = cell2mat(DF.Po);
DFVB = cell2mat(DF.VB);
DFPoVPM = cell2mat(DF.PoVPM);
DFLD = cell2mat(DF.LD);
DFnRT = cell2mat(DF.nRT);
[Wp(1) Wh(1)] = b_ranksum2(-DFPo,zeros(size(DFPo)),'alpha',0.05);   % Mann-Whitney test for H0: all diff >= 0
                                                                    % against H1: all diff < 0
[Wp(2) Wh(2)] = b_ranksum2(-DFVB,zeros(size(DFVB)),'alpha',0.05);
[Wp(3) Wh(3)] = b_ranksum2(-DFPoVPM,zeros(size(DFPoVPM)),'alpha',0.05);
[Wp(4) Wh(4)] = b_ranksum2(-DFLD,zeros(size(DFLD)),'alpha',0.05);
[Wp(5) Wh(5)] = b_ranksum2(-DFnRT,zeros(size(DFnRT)),'alpha',0.05);
dfstat.Po.allno = length(DFPo);     % number of differnces
dfstat.VB.allno = length(DFVB);
dfstat.PoVPM.allno = length(DFPoVPM);
dfstat.LD.allno = length(DFLD);
dfstat.nRT.allno = length(DFnRT);
dfstat.Po.posno = length(find(DFPo>0));     % number of pos. differnces
dfstat.VB.posno = length(find(DFVB>0));
dfstat.PoVPM.posno = length(find(DFPoVPM>0));
dfstat.LD.posno = length(find(DFLD>0));
dfstat.nRT.posno = length(find(DFnRT>0));
dfstat.Po.posprop = dfstat.Po.posno / dfstat.Po.allno;     % proportion of pos. differnces
dfstat.VB.posprop = dfstat.VB.posno / dfstat.VB.allno;
dfstat.PoVPM.posprop = dfstat.PoVPM.posno / dfstat.PoVPM.allno;
dfstat.LD.posprop = dfstat.LD.posno / dfstat.LD.allno;
dfstat.nRT.posprop = dfstat.nRT.posno / dfstat.nRT.allno;
dfstat.Po.dfmean = circular_mean(DFPo,'deg');       % overall mean difference
dfstat.VB.dfmean = circular_mean(DFVB,'deg');
dfstat.PoVPM.dfmean = circular_mean(DFPoVPM,'deg');
dfstat.LD.dfmean = circular_mean(DFLD,'deg');
dfstat.nRT.dfmean = circular_mean(DFnRT,'deg');
dfstat.Po.alldiffp = Wp(1);   % Mann-Whitney test for H0: all diff >= 0 against H1: all diff < 0
dfstat.VB.alldiffp = Wp(2);
dfstat.PoVPM.alldiffp = Wp(3);
dfstat.LD.alldiffp = Wp(4);
dfstat.nRT.alldiffp = Wp(5);
cd(resdir)
save dfstat dfstat

% Plot & save
dbclear if error
cd(resdir)
H = figure;
boxplot([DFPo DFVB DFPoVPM DFLD DFnRT],[zeros(size(DFPo)) ones(size(DFVB))...
    2*ones(size(DFPoVPM)) 3*ones(size(DFLD)) 4*ones(size(DFnRT))],'Labels',...
    {'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
y_lim = ylim;
x_lim = xlim;
for k = 1:5
    tpos1 = k;
    tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
    if Wh(k)
        clr = 'red';
    else
        clr = 'black';
    end
    text(tpos1,tpos2,num2str(Wp(k)),'Color',clr,'Horizontalalignment','center')
end
line([x_lim(1) x_lim(2)],[0 0],'Color','green')
title('Burstphase differences pooled')
saveas(H,'burstphase6_box.fig')
saveas(H,'burstphase6_box.eps')

DFpooled = {DFPo DFVB DFPoVPM DFLD DFnRT};
H = dfplot(DFpooled);
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
title('Burstphase MRL differences pooled with SE')
saveas(H,'burstphase7_pooled.fig')
saveas(H,'burstphase7_pooled.eps')

H = dfplot2(DFpooled);
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
title('Burstphase MRL differences pooled with all values')
saveas(H,'burstphase8_pooled.fig')
saveas(H,'burstphase8_pooled.eps')

H = dfplot3(DFpooled);
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
title('Burstphase MRL differences pooled with conf. int.')
saveas(H,'burstphase9_pooled.fig')
saveas(H,'burstphase9_pooled.eps')

H = dfplot(DF.Po);
title('Po burstphase differences with SE')
saveas(H,'burstphase3_Po.fig')
saveas(H,'burstphase3_Po.eps')
H = dfplot2(DF.Po);
title('Po burstphase differences with all values')
saveas(H,'burstphase4_Po.fig')
saveas(H,'burstphase4_Po.eps')
H = dfplot3(DF.Po);
title('Po burstphase differences with conf. int.')
saveas(H,'burstphase5_Po.fig')
saveas(H,'burstphase5_Po.eps')

H = dfplot(DF.VB);
title('VB burstphase differences with SE')
saveas(H,'burstphase3_VB.fig')
saveas(H,'burstphase3_VB.eps')
H = dfplot2(DF.VB);
title('VB burstphase differences with all values')
saveas(H,'burstphase4_VB.fig')
saveas(H,'burstphase4_VB.eps')
H = dfplot3(DF.VB);
title('VB burstphase differences with conf. int.')
saveas(H,'burstphase5_VB.fig')
saveas(H,'burstphase5_VB.eps')

H = dfplot(DF.PoVPM);
title('PoVPM burstphase differences with SE')
saveas(H,'burstphase3_PoVPM.fig')
saveas(H,'burstphase3_PoVPM.eps')
H = dfplot2(DF.PoVPM);
title('PoVPM burstphase differences with all values')
saveas(H,'burstphase4_PoVPM.fig')
saveas(H,'burstphase4_PoVPM.eps')
H = dfplot3(DF.PoVPM);
title('PoVPM burstphase differences with conf. int.')
saveas(H,'burstphase5_PoVPM.fig')
saveas(H,'burstphase5_PoVPM.eps')

H = dfplot(DF.LD);
title('LD burstphase differences with SE')
saveas(H,'burstphase3_LD.fig')
saveas(H,'burstphase3_LD.eps')
H = dfplot2(DF.LD);
title('LD burstphase differences with all values')
saveas(H,'burstphase4_LD.fig')
saveas(H,'burstphase4_LD.eps')
H = dfplot3(DF.LD);
title('LD burstphase differences with conf. int.')
saveas(H,'burstphase5_LD.fig')
saveas(H,'burstphase5_LD.eps')

H = dfplot(DF.nRT);
title('nRT burstphase differences with SE')
saveas(H,'burstphase3_nRT.fig')
saveas(H,'burstphase3_nRT.eps')
H = dfplot2(DF.nRT);
title('nRT burstphase differences with all values')
saveas(H,'burstphase4_nRT.fig')
saveas(H,'burstphase4_nRT.eps')
H = dfplot3(DF.nRT);
title('nRT burstphase differences with conf. int.')
saveas(H,'burstphase5_nRT.fig')
saveas(H,'burstphase5_nRT.eps')

close all
cd(mm)

% -------------------------------------------------------------------------
function [ftm, ang, mvl] = mvlmn(ng,degrad)

if isequal(degrad,'deg')    % convert to radian
    ng = ng / 180 * pi;
end
ftm = sum(exp(1).^(i*ng)) / length(ng);    % first trigonometric moment
ang = angle(ftm);   % mean angle
mvl = abs(ftm);     % mean resultant length
if isequal(degrad,'deg')    % convert to degrees
    ang = ang * 180 / pi;
end

% -------------------------------------------------------------------------
function H = dfplot(df)

ldf = length(df);
mndf = zeros(1,ldf);
SE = zeros(1,ldf);
for k = 1:ldf
    [ftm, ang, mvl] = mvlmn(df{k},'deg');
    mndf(k) = ang;
    SE(k) = circular_SE(df{k},'deg');
end
H = figure;
bar((1:ldf),mndf)
hold on
errorbar((1:ldf),mndf,SE,'k+')

% -------------------------------------------------------------------------
function H = dfplot2(df)

ldf = length(df);
mndf = zeros(1,ldf);
for k = 1:ldf
    [ftm, ang, mvl] = mvlmn(df{k},'deg');
    mndf(k) = ang;
end
H = figure;
bar((1:ldf),mndf,'FaceColor','none')
hold on
for k = 1:ldf
    plot(k*ones(size(df{k})),df{k},'.')
end

% -------------------------------------------------------------------------
function H = dfplot3(df)

ldf = length(df);
mndf = zeros(1,ldf);
L = zeros(1,ldf);
U = zeros(1,ldf);
for k = 1:ldf
    [ftm, ang, mvl] = mvlmn(df{k},'deg');
    mndf(k) = ang;
    [L(k) U(k)] = circular_conf(df{k},'deg');
end
H = figure;
bar((1:ldf),mndf)
hold on
errorbar((1:ldf),mndf,mndf-L,'k+')

% -------------------------------------------------------------------------
function df = circdiff(P2,P1)
% P2 - P1 on the circle

df1 = P2 - P1;
if df1 < 0
    df2 = 360 + df1;
else
    df2 = df1 - 360;
end
mdf = min(abs(df1),abs(df2));
if abs(df1) == mdf
    df = df1;
elseif abs(df2) == mdf
    df = df2;
elseif isnan(mdf)
    df = NaN;
else
    error('Technical error 197')
end