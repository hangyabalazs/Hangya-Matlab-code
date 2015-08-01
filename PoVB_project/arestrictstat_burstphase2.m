function arestrictstat_burstphase2
%ARESTRICTSTAT_BURSTPHASE2   Summary phase statistics of ordinal intraburst spikes in bounded frequency band.
%   ARESTRICTSTAT_BURSTPHASE2 calculates summary phase statistics of
%   ordinal intraburst spikes for Po, VPM, PoVPM, LD and nRT data. It plots
%   and saves result figures displaying mean ordinal spike phases with
%   circular standard error bars. Edit code to modify input and output
%   directories!
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
ANG = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
HPo = figure;
title('Po')
hold on
HVB = figure;
title('VB')
hold on
HPoVPM = figure;
title('PoVPM')
hold on
HLD = figure;
title('LD')
hold on
HnRT = figure;
title('nRT')
hold on
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
    SE = zeros(1,mibs-1);
    aang_ord_all = cell(1,mibs-1);
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
            SE(k-1) = circular_SE(ng,'rad');
            nt = eval(['length(MRL.' loc ');']);
            if k - 1 > nt
                ntt = 0;
                eval(['ANG.' loc '{k-1} = [];']);
            else
                ntt = eval(['length(MRL.' loc '{k-1});']);
            end
            eval(['MRL.' loc '{k-1}(ntt+1) = mrl(k-1);']);
            eval(['MN.' loc '{k-1}(ntt+1) = mn(k-1);']);
            aang_ord_all{k-1} = ng;
        else
            mn_rad(k-1) = NaN;
            mn(k-1) = NaN;
            mrl(k-1) = NaN;
            SE(k-1) = NaN;
        end
    end
    
    cd(resdir)      % plot
    dbclear if error
    H = polplot(aang_ord_all);  
    title(gca,[tstr ' ' loc])
    saveas(H,[fname '_' loc])
    close(H)
    
    figstr = ['H' loc];
    figure(eval(figstr))
    linplot(mn,SE)
    
    dbstop if error
end

% Save
dbclear if error
cd(resdir)
saveas(HPo,'burstphase_Po.fig')
saveas(HPo,'burstphase_Po.eps')
saveas(HVB,'burstphase_VB.fig')
saveas(HVB,'burstphase_VB.eps')
saveas(HPoVPM,'burstphase_PoVPM.fig')
saveas(HPoVPM,'burstphase_PoVPM.eps')
saveas(HLD,'burstphase_LD.fig')
saveas(HLD,'burstphase_LD.eps')
saveas(HnRT,'burstphase_nRT.fig')
saveas(HnRT,'burstphase_nRT.eps')

cd(mm)

% -------------------------------------------------------------------------
function ftm = mvlmn(fr1_all)

ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment

% -------------------------------------------------------------------------
function H = polplot(aa)

clr = {'red' 'green' 'blue' 'magenta' 'yellow' 'cyan' ...
        [0 153 0]/256 [255 102 0]/256 [255 204 204]/256 [102 102 102]/256};
H = figure;
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white','LineWidth',0.01)
hold on
CMP = [];
legstr = {};
for k = 1:length(aa)
    if ~isempty(aa{k})
        lftm = mvlmn(aa{k});
    else
        lftm = 0;
    end
    CMP(end+1) = compass(real(lftm),imag(lftm));
    set(CMP(end),'Color',clr{k},'LineWidth',2)
    legstr{end+1} = num2str(k+1);
end
[L,OBJ] = legend(CMP,legstr,'Location','NorthEastOutside');

% -------------------------------------------------------------------------
function linplot(mn,SE)

mibs = length(mn) + 1;
plot([2:mibs],mn)
for k = 2:mibs
    if ~isnan(mn(k-1))
        line([k k],[mn(k-1)-SE(k-1)/pi*180 mn(k-1)+SE(k-1)/pi*180],'Color','red')
    end
end
x_lim = xlim;
xl2 = max(x_lim(2),mibs);
xlim([1 xl2])