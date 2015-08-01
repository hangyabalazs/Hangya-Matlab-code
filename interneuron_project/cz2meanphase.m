function cz2meanphase
%CZ2MEANPHASE   Mean phase.
%   CZ2MEANPHASE calculates mean wavelet phase and mean phase of the
%   original calculations by Yu Li. Results are saved in an Excel table.
%   Note that phase values are corrected in a way that zero phase
%   corresponds to the negative peaks of the theta cycle!
%
%   See also CZ2PHASE.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Czurko2\Phase2\wavelet2\angs\'];
table = [DATAPATH 'Czurko2\orig_phases.xls'];
resdir = [DATAPATH 'Czurko2\Phase2\'];
files = b_filelist(inpdir);
sf = length(files);
mm = pwd;
cd(resdir)

% Progress indicator
wb = waitbar(0,'Running CZ2MEANPHASE...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Load EEG table
[tbl0 tbl00 tbl] = xlsread(table,'Sheet1');

% Main: calculate phase distribution
xlsout = cell(sf,8);    % initialize Excel output
for o = 2:sf
    disp(o)
    fname = files(o).name
    fs = findstr(fname,'_');
    fn = fname(1:fs(1)-1);
    fnp = fname(fs(1)+1:fs(2)-1);
    fnc = fname(fs(2)+1:fs(3)-1);
    ff = [inpdir fname];
    load(ff)   % load wavelet phase values
    
    wdeg = wang * 180 / pi;   % mean of wavelet phase
    wdeg = wdeg - 180;        % zero phase: negative peak
    wdeg(wdeg<-180) = wdeg(wdeg<-180) + 360;
    
    inx = find(strcmp(cellfun(@eegname,{tbl{1:end,2}},'UniformOutput',0),fn));
    if length(inx) > 1
        inx2 = find(strcmp({tbl{inx,3}},upper(fnp))&strcmp({tbl{inx,4}},upper(fnc)));
        if length(inx2) > 1
            error('Technical error 49.')
        end
        inx = inx(inx2);
    end
        
    orig_phs = tbl{inx,6};  % load original phase values
    phs2 = strread(orig_phs,'%s','delimiter',' ');
    phs3 = cellfun(@str2num,phs2)';
    orig_maxphase = tbl{inx,7};
    
    bs = [280:20:340 0:20:260];   % mean of original phase
    cuvs = exp(i*(bs/180*pi));
    wcuvs = cuvs .* phs3;
    mv = sum(wcuvs);
    mang = angle(mv);
    mangdeg = mang * 180 / pi;
    mangdeg(mangdeg<0) = mangdeg(mangdeg<0) + 360;
    
%     bs2 = [0:20:340];
%     edges = [10:20:330];
    bs2 = [0:18:342];
    wangs2 = wangs * 180 / pi;
    wangs2(wangs2<-10) = wangs2(wangs2<-10) + 360;
    xout = hist(wangs2,bs2);      % maximal firing phase (wavelet)
%     xout = histc(wangs2,edges);
%     xout = xout(1:end-1);
    xxout = [sum(xout(1:3)) sum(xout(2:4)) sum(xout(3:5)) sum(xout(4:6)) sum(xout(5:7)) sum(xout(6:8)) sum(xout(7:9)) sum(xout(8:10)) ...
        sum(xout(9:11)) sum(xout(10:12)) sum(xout(11:13)) sum(xout(12:14)) sum(xout(13:15)) sum(xout(14:16)) ...
        sum(xout(15:17)) sum(xout(16:18)) sum([xout(17:18) xout(1)]) sum([xout(18) xout(1:2)])];
    fx = find(xxout==max(xxout)) + 1;
    fx = mod2(fx,18);
    wmax = bs2(fx) - 180;
    wmax = circular_mean(wmax,'deg');
    wmax(wmax<0) = wmax(wmax<0) + 360;
    
%     xxout = [sum(phs3(1:3)) sum(phs3(2:4)) sum(phs3(3:5)) sum(phs3(4:6)) sum(phs3(5:7)) sum(phs3(6:8)) sum(phs3(7:9)) sum(phs3(8:10)) ...
%         sum(phs3(9:11)) sum(phs3(10:12)) sum(phs3(11:13)) sum(phs3(12:14)) sum(phs3(13:15)) sum(phs3(14:16)) ...
%         sum(phs3(15:17)) sum(phs3(16:18)) sum([phs3(17:18) phs3(1)]) sum([phs3(18) phs3(1:2)])];
%     fx = find(xxout==max(xxout)) + 1;   % mean of three largest bin (original)
%     fx = mod2(fx,18);
    fx = find(bs(1:4)<=orig_maxphase,1,'last'); 
    if isempty(fx)
        fx = 4 + find(bs(5:end)<=orig_maxphase,1,'last');
    end
    if mod2(bs(mod2(fx+1,18)),360) - orig_maxphase < orig_maxphase - bs(fx)
        fx = fx + 1;
        fx = mod2(fx,18);
    end
    cuvs = exp(i*([bs(mod2(fx-1,18)) bs(fx) bs(mod2(fx+1,18))]/180*pi));
    wcuvs = cuvs .* [phs3(mod2(fx-1,18)) phs3(fx) phs3(mod2(fx+1,18))];
    mv = sum(wcuvs);
    mang = angle(mv);
    mangdeg2 = mang * 180 / pi;
    mangdeg2(mangdeg2<0) = mangdeg2(mangdeg2<0) + 360;
    omax = bs(fx);   % maximal firing phase (original)
    omax(omax<0) = omax(omax<0) + 360;
    
    t = 2;
    cuvs = exp(i*(bs(mod2(fx-t:fx+t,18))/180*pi));
    wcuvs = cuvs .* phs3(mod2(fx-t:fx+t,18));
    mv = sum(wcuvs);
    mang = angle(mv);
    mangdeg3 = mang * 180 / pi;
    mangdeg3(mangdeg3<0) = mangdeg3(mangdeg3<0) + 360;
        
    xlsout{o,1} = fname;    % Excel output
    xlsout{o,2} = wdeg;
    xlsout{o,3} = mangdeg;
    xlsout{o,4} = wmax;
    xlsout{o,5} = mangdeg2;
    xlsout{o,6} = omax;
    xlsout{o,7} = orig_maxphase;
    xlsout{o,8} = mangdeg3;
end
xlswrite('mean_phases2',xlsout)   % write Excel output
close(wb)
cd(mm)

% -------------------------------------------------------------------------
function A = eegname(S)

fs = findstr(S,'/');
if isempty(fs)
    A = '';
    return
end
A = S(fs(end)+1:end);