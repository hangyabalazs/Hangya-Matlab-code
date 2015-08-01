function arestrictstat_burstphase4
%ARESTRICTSTAT_BURSTPHASE4   Summary phase statistics of ordinal intraburst spikes in bounded frequency band.
%   ARESTRICTSTAT_BURSTPHASE4 calculates summary phase statistics of
%   ordinal intraburst spikes for Po, VPM, PoVPM, LD and nRT data. It plots
%   and saves result figures displaying mean resultant length of ordinal
%   spike phases. Edit code to modify input and output directories!
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
    end
end

% Plot & save
cd(resdir)
leg_handle1 = [];
leg_handle2 = [];
H1 = figure;
hold on
H2 = figure;
hold on
[leg_handle1(end+1) leg_handle2(end+1)] = mrlplot(H1,H2,MRL.Po,'b','o');
[leg_handle1(end+1) leg_handle2(end+1)] = mrlplot(H1,H2,MRL.VB,'r','v');
[leg_handle1(end+1) leg_handle2(end+1)] = mrlplot(H1,H2,MRL.PoVPM,'m','+');
[leg_handle1(end+1) leg_handle2(end+1)] = mrlplot(H1,H2,MRL.LD,'g','d');
[leg_handle1(end+1) leg_handle2(end+1)] = mrlplot(H1,H2,MRL.nRT,'y','x');
figure(H1)
legend(leg_handle1,{'Po' 'VB' 'PoVPM' 'LD' 'nRT'})
title('Mean+SE MRL of ordinal spkike distributions')
figure(H2)
legend(leg_handle2,{'Po' 'VB' 'PoVPM' 'LD' 'nRT'})
title('Mean+SD MRL of ordinal spkike distributions')
dbclear if error
saveas(H1,'burstphase_meanmrl1.fig')
saveas(H1,'burstphase_meanmrl1.eps')
saveas(H2,'burstphase_meanmrl2.fig')
saveas(H2,'burstphase_meanmrl2.eps')

cd(resdir)
H = Rplot(MRL.Po);
title('Po burstphase MVLs with SE')
saveas(H,'burstphase_meanmrl3_Po.fig')
saveas(H,'burstphase_meanmrl3_Po.eps')
H = Rplot2(MRL.Po);
title('Po burstphase MVLs with all values')
saveas(H,'burstphase_meanmrl4_Po.fig')
saveas(H,'burstphase_meanmrl4_Po.eps')
H = Rplot3(MRL.Po);
title('Po burstphase MVLs with SD')
saveas(H,'burstphase_meanmrl5_Po.fig')
saveas(H,'burstphase_meanmrl5_Po.eps')
H = Rbox(MRL.Po);
title('Po burstphase MVL boxplots')
saveas(H,'burstphase_meanmrl6_Po.fig')
saveas(H,'burstphase_meanmrl6_Po.eps')

H = Rplot(MRL.VB);
title('VB burstphase MVSs with SE')
saveas(H,'burstphase_meanmrl3_VB.fig')
saveas(H,'burstphase_meanmrl3_VB.eps')
H = Rplot2(MRL.VB);
title('VB burstphase MVLs with all values')
saveas(H,'burstphase_meanmrl4_VB.fig')
saveas(H,'burstphase_meanmrl4_VB.eps')
H = Rplot3(MRL.VB);
title('VB burstphase MVLs with SD')
saveas(H,'burstphase_meanmrl5_VB.fig')
saveas(H,'burstphase_meanmrl5_VB.eps')
H = Rbox(MRL.VB);
title('VB burstphase MVL boxplots')
saveas(H,'burstphase_meanmrl6_VB.fig')
saveas(H,'burstphase_meanmrl6_VB.eps')

H = Rplot(MRL.PoVPM);
title('PoVPM burstphase MVLs with SE')
saveas(H,'burstphase_meanmrl3_PoVPM.fig')
saveas(H,'burstphase_meanmrl3_PoVPM.eps')
H = Rplot2(MRL.PoVPM);
title('PoVPM burstphase MVLs with all values')
saveas(H,'burstphase_meanmrl4_PoVPM.fig')
saveas(H,'burstphase_meanmrl4_PoVPM.eps')
H = Rplot3(MRL.PoVPM);
title('PoVPM burstphase MVLs with SD')
saveas(H,'burstphase_meanmrl5_PoVPM.fig')
saveas(H,'burstphase_meanmrl5_PoVPM.eps')
H = Rbox(MRL.PoVPM);
title('PoVPM burstphase MVL boxplots')
saveas(H,'burstphase_meanmrl6_PoVPM.fig')
saveas(H,'burstphase_meanmrl6_PoVPM.eps')

H = Rplot(MRL.LD);
title('LD burstphase MVLs with SE')
saveas(H,'burstphase_meanmrl3_LD.fig')
saveas(H,'burstphase_meanmrl3_LD.eps')
H = Rplot2(MRL.LD);
title('LD burstphase MVLs with all values')
saveas(H,'burstphase_meanmrl4_LD.fig')
saveas(H,'burstphase_meanmrl4_LD.eps')
H = Rplot3(MRL.LD);
title('LD burstphase MVLs with SD')
saveas(H,'burstphase_meanmrl5_LD.fig')
saveas(H,'burstphase_meanmrl5_LD.eps')
H = Rbox(MRL.LD);
title('LD burstphase MVL boxplots')
saveas(H,'burstphase_meanmrl6_LD.fig')
saveas(H,'burstphase_meanmrl6_LD.eps')

H = Rplot(MRL.nRT);
title('nRT burstphase MVLs with SE')
saveas(H,'burstphase_meanmrl3_nRT.fig')
saveas(H,'burstphase_meanmrl3_nRT.eps')
H = Rplot2(MRL.nRT);
title('nRT burstphase MVLs with all values')
saveas(H,'burstphase_meanmrl4_nRT.fig')
saveas(H,'burstphase_meanmrl4_nRT.eps')
H = Rplot3(MRL.nRT);
title('nRT burstphase MVLs with SD')
saveas(H,'burstphase_meanmrl5_nRT.fig')
saveas(H,'burstphase_meanmrl5_nRT.eps')
H = Rbox(MRL.nRT);
title('nRT burstphase MVL boxplots')
saveas(H,'burstphase_meanmrl6_nRT.fig')
saveas(H,'burstphase_meanmrl6_nRT.eps')
% close all
cd(mm)

% -------------------------------------------------------------------------
function [L1, L2] = mrlplot(H1,H2,mrl,color,mark)

lmrl = length(mrl);
SE = zeros(1,lmrl);
SD = zeros(1,lmrl);
mean_mrl = zeros(1,lmrl);
for k = 1:lmrl
    lm = length(mrl{k});
    SE(k) = std(mrl{k}) / sqrt(lm);
    SD(k) = std(mrl{k});
    mean_mrl(k) = mean(mrl{k});
end
lsm1 = [color '-' mark];      % line style, marker
lsm2 = [color '+'];
figure(H1)
L1 = plot(2:lmrl+1,mean_mrl,lsm1);   % handle for the legend
hold on
errorbar((2:lmrl+1),mean_mrl,SE,lsm2)
figure(H2)
L2 = plot(2:lmrl+1,mean_mrl,lsm1);
hold on
errorbar((2:lmrl+1),mean_mrl,SD,lsm2)

% -------------------------------------------------------------------------
function H = Rplot(mvl)

lmvl = length(mvl);
mnmvl = zeros(1,lmvl);
SE = zeros(1,lmvl);
for k = 1:lmvl
    mnmvl(k) = mean(mvl{k});
    SE(k) = std(mvl{k}) / sqrt(length(mvl{k}));
end
H = figure;
bar((2:lmvl+1),mnmvl)
hold on
errorbar((2:lmvl+1),mnmvl,SE,'k+')

% -------------------------------------------------------------------------
function H = Rplot2(mvl)

lmvl = length(mvl);
mnmvl = zeros(1,lmvl);
for k = 1:lmvl
    mnmvl(k) = mean(mvl{k});
end
H = figure;
bar((2:lmvl+1),mnmvl,'FaceColor','none')
hold on
for k = 1:lmvl
    plot((k+1)*ones(size(mvl{k})),mvl{k},'.')
end

% -------------------------------------------------------------------------
function H = Rplot3(mvl)

lmvl = length(mvl);
mnmvl = zeros(1,lmvl);
SD = zeros(1,lmvl);
for k = 1:lmvl
    mnmvl(k) = mean(mvl{k});
    SD(k) = std(mvl{k});
end
H = figure;
bar((2:lmvl+1),mnmvl)
hold on
errorbar((2:lmvl+1),mnmvl,SD,'k+')

% -------------------------------------------------------------------------
function H = Rbox(mvl)

lmvl = length(mvl);
m = [];
g = [];     % grouping variable
l = {};     % labels
Wp = [];    % Mann-Whitney test
clr = {};
for k = 1:lmvl
    mt = mvl{k};
    m = [m mt];
    g = [g (k-1)*ones(size(mt))];
    l = [l num2str(k+1)];
    if k < lmvl
        mtt = mvl{k+1};
        [Wp(end+1),Wh] = b_ranksum2(mt,mtt,'alpha',0.05);
        if Wh
            clr{end+1} = 'red';
        else
            clr{end+1} = 'black';
        end
    end
end
H = figure;
boxplot(m,g,'labels',l);

y_lim = ylim;
x_lim = xlim;
for k = 1:lmvl-1
    tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * k / lmvl;
    tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
    text(tpos1,tpos2,num2str(Wp(k)),'Color',clr{k},'Horizontalalignment','center')
end