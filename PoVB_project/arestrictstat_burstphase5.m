function arestrictstat_burstphase5
%ARESTRICTSTAT_BURSTPHASE5   Summary phase statistics of ordinal intraburst spikes in bounded frequency band.
%   ARESTRICTSTAT_BURSTPHASE5 calculates summary phase statistics of
%   ordinal intraburst spikes for Po, VPM, PoVPM, LD and nRT data. It plots
%   and saves result figures displaying MVL differences of ordinal spike
%   phases. A structure ('dfmvlstat') containing the number and proportion
%   of negative differences as well as the overall mean and median
%   difference is also saved. Boxplot and bargraphs of pooled differences
%   are saved and their relation to zero median is tested using
%   Mann-Whitney U-test. Test results appear on the boxplot as well as in
%   the 'dfmvlstat' structure. Edit code to modify input and output
%   directories!
%
%   See also ARESTRICTSTAT_BURST, ARESTRICTSTAT_PHASE3 and
%   ARESTRICTSTAT_BURSTPHASE3.

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
DFMRL = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
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
    dfmrl = zeros(1,mibs-2);
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
            dfmrl(k-2) = mrl(k-1) - mrl(k-2);
        end
    end
    for k = 2:mibs-1
        nt = eval(['length(DFMRL.' loc ');']);
        if k - 1 > nt
            ntt = 0;
        else
            ntt = eval(['length(DFMRL.' loc '{k-1});']);
        end
        if ~isnan(dfmrl(k-1))
            eval(['DFMRL.' loc '{k-1}(ntt+1) = dfmrl(k-1);']);
        end
    end
end

% Saving statistics
DFPo = cell2mat(DFMRL.Po);
DFVB = cell2mat(DFMRL.VB);
DFPoVPM = cell2mat(DFMRL.PoVPM);
DFLD = cell2mat(DFMRL.LD);
DFnRT = cell2mat(DFMRL.nRT);
[Wp(1) Wh(1)] = b_ranksum2(DFPo,zeros(size(DFPo)),'alpha',0.05);   % Mann-Whitney test for H0: all diff <= 0
                                                                    % against H1: all diff > 0
[Wp(2) Wh(2)] = b_ranksum2(DFVB,zeros(size(DFVB)),'alpha',0.05);
[Wp(3) Wh(3)] = b_ranksum2(DFPoVPM,zeros(size(DFPoVPM)),'alpha',0.05);
[Wp(4) Wh(4)] = b_ranksum2(DFLD,zeros(size(DFLD)),'alpha',0.05);
[Wp(5) Wh(5)] = b_ranksum2(DFnRT,zeros(size(DFnRT)),'alpha',0.05);
dfmvlstat.Po.allno = length(DFPo);     % number of differnces
dfmvlstat.VB.allno = length(DFVB);
dfmvlstat.PoVPM.allno = length(DFPoVPM);
dfmvlstat.LD.allno = length(DFLD);
dfmvlstat.nRT.allno = length(DFnRT);
dfmvlstat.Po.negno = length(find(DFPo<0));     % number of neg. differnces
dfmvlstat.VB.negno = length(find(DFVB<0));
dfmvlstat.PoVPM.negno = length(find(DFPoVPM<0));
dfmvlstat.LD.negno = length(find(DFLD<0));
dfmvlstat.nRT.negno = length(find(DFnRT<0));
dfmvlstat.Po.negprop = dfmvlstat.Po.negno / dfmvlstat.Po.allno;     % proportion of neg. differnces
dfmvlstat.VB.negprop = dfmvlstat.VB.negno / dfmvlstat.VB.allno;
dfmvlstat.PoVPM.negprop = dfmvlstat.PoVPM.negno / dfmvlstat.PoVPM.allno;
dfmvlstat.LD.negprop = dfmvlstat.LD.negno / dfmvlstat.LD.allno;
dfmvlstat.nRT.negprop = dfmvlstat.nRT.negno / dfmvlstat.nRT.allno;
dfmvlstat.Po.dfmvlmean = mean(DFPo);       % overall mean difference
dfmvlstat.VB.dfmvlmean = mean(DFVB);
dfmvlstat.PoVPM.dfmvlmean = mean(DFPoVPM);
dfmvlstat.LD.dfmvlmean = mean(DFLD);
dfmvlstat.nRT.dfmvlmean = mean(DFnRT);
dfmvlstat.Po.dfmvlmedian = median(DFPo);       % overall median difference
dfmvlstat.VB.dfmvlmedian = median(DFVB);
dfmvlstat.PoVPM.dfmvlmedian = median(DFPoVPM);
dfmvlstat.LD.dfmvlmedian = median(DFLD);
dfmvlstat.nRT.dfmvlmedian = median(DFnRT);
dfmvlstat.Po.alldiffp = Wp(1);   % Mann-Whitney test for H0: all diff <= 0 against H1: all diff > 0
dfmvlstat.VB.alldiffp = Wp(2);
dfmvlstat.PoVPM.alldiffp = Wp(3);
dfmvlstat.LD.alldiffp = Wp(4);
dfmvlstat.nRT.alldiffp = Wp(5);
cd(resdir)
save dfmvlstat dfmvlstat

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
title('Burstphase MRL differences pooled')
saveas(H,'burstphase_dfmrl1_box.fig')
saveas(H,'burstphase_dfmrl1_box.eps')

DFpooled = {DFPo DFVB DFPoVPM DFLD DFnRT};
H = dfplot(DFpooled);
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
title('Burstphase MRL differences pooled with SE')
saveas(H,'burstphase_dfmrl5.fig')
saveas(H,'burstphase_dfmrl5.eps')

H = dfplot2(DFpooled);
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
title('Burstphase MRL differences pooled with all values')
saveas(H,'burstphase_dfmrl6.fig')
saveas(H,'burstphase_dfmrl6.eps')

H = dfplot3(DFpooled);
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'Po' 'VB' 'PoVPM' 'LD' 'nRT'});
title('Burstphase MRL differences pooled with SD')
saveas(H,'burstphase_dfmrl7.fig')
saveas(H,'burstphase_dfmrl7.eps')

H = dfplot(DFMRL.Po);
title('Po burstphase MRL differences with SE')
saveas(H,'burstphase_dfmrl2_Po.fig')
saveas(H,'burstphase_dfmrl2_Po.eps')
H = dfplot2(DFMRL.Po);
title('Po burstphase MRL differences with all values')
saveas(H,'burstphase_dfmrl3_Po.fig')
saveas(H,'burstphase_dfmrl3_Po.eps')
H = dfplot3(DFMRL.Po);
title('Po burstphase MRL differences with SD')
saveas(H,'burstphase_dfmrl4_Po.fig')
saveas(H,'burstphase_dfmrl4_Po.eps')

H = dfplot(DFMRL.VB);
title('VB burstphase MRL differences with SE')
saveas(H,'burstphase_dfmrl2_VB.fig')
saveas(H,'burstphase_dfmrl2_VB.eps')
H = dfplot2(DFMRL.VB);
title('VB burstphase MRL differences with all values')
saveas(H,'burstphase_dfmrl3_VB.fig')
saveas(H,'burstphase_dfmrl3_VB.eps')
H = dfplot3(DFMRL.VB);
title('VB burstphase MRL differences with SD')
saveas(H,'burstphase_dfmrl4_VB.fig')
saveas(H,'burstphase_dfmrl4_VB.eps')

H = dfplot(DFMRL.PoVPM);
title('PoVPM burstphase MRL differences with SE')
saveas(H,'burstphase_dfmrl2_PoVPM.fig')
saveas(H,'burstphase_dfmrl2_PoVPM.eps')
H = dfplot2(DFMRL.PoVPM);
title('PoVPM burstphase MRL differences with all values')
saveas(H,'burstphase_dfmrl3_PoVPM.fig')
saveas(H,'burstphase_dfmrl3_PoVPM.eps')
H = dfplot3(DFMRL.PoVPM);
title('PoVPM burstphase MRL differences with SD')
saveas(H,'burstphase_dfmrl4_PoVPM.fig')
saveas(H,'burstphase_dfmrl4_PoVPM.eps')

H = dfplot(DFMRL.LD);
title('LD burstphase MRL differences with SE')
saveas(H,'burstphase_dfmrl2_LD.fig')
saveas(H,'burstphase_dfmrl2_LD.eps')
H = dfplot2(DFMRL.LD);
title('LD burstphase MRL differences with all values')
saveas(H,'burstphase_dfmrl3_LD.fig')
saveas(H,'burstphase_dfmrl3_LD.eps')
H = dfplot3(DFMRL.LD);
title('LD burstphase MRL differences with SD')
saveas(H,'burstphase_dfmrl4_LD.fig')
saveas(H,'burstphase_dfmrl4_LD.eps')

H = dfplot(DFMRL.nRT);
title('nRT burstphase MRL differences with SE')
saveas(H,'burstphase_dfmrl2_nRT.fig')
saveas(H,'burstphase_dfmrl2_nRT.eps')
H = dfplot2(DFMRL.nRT);
title('nRT burstphase MRL differences with all values')
saveas(H,'burstphase_dfmrl3_nRT.fig')
saveas(H,'burstphase_dfmrl3_nRT.eps')
H = dfplot3(DFMRL.nRT);
title('nRT burstphase MRL differences with SD')
saveas(H,'burstphase_dfmrl4_nRT.fig')
saveas(H,'burstphase_dfmrl4_nRT.eps')

% close all
cd(mm)

% -------------------------------------------------------------------------
function H = dfplot(df)

ldf = length(df);
mndf = zeros(1,ldf);
SE = zeros(1,ldf);
for k = 1:ldf
    mndf(k) = mean(df{k});
    SE(k) = std(df{k}) / sqrt(length(df{k}));
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
    mndf(k) = mean(df{k});
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
SD = zeros(1,ldf);
for k = 1:ldf
    mndf(k) = mean(df{k});
    SD(k) = std(df{k});
end
H = figure;
bar((1:ldf),mndf)
hold on
errorbar((1:ldf),mndf,SD,'k+')