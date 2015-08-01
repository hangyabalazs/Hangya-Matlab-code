function rapheconfocal_run
%RAPHECONFOCAL_RUN   Pixel counts on a sequence of files.
%   RAPHECONFOCAL_RUN calculates pixel counts on a sequence of files (see
%   RAPHECONFOCAL for details). Pixel counts are averaged over the images
%   (average is triggered by the center of the pyramidal layer). The ratio
%   of overlapping points to green pixels is also calculated, plotted and
%   saved. Lines are smoothed using a sliding average with 20 um windowsize.
%   Thresholded (binary) images are saved.
%
%   See also RAPHECONFOCAL, RAPHEPYRLAY and RAPHERELPYRPOS.

% Directories
global DATAPATH
cnm = '875';
inpdir_gfp = ['X:\Zsolt\zs7\zs7g_rost_denzitas\balazs\' cnm '\gfp\'];
inpdir_ser = ['X:\Zsolt\zs7\zs7g_rost_denzitas\balazs\' cnm '\ser\'];
inpdir_pyrpos = [DATAPATH 'Raphe\zs7\pyrlay\'];
inpdir_relpyrpos = [DATAPATH 'Raphe\zs7\relpyrpos\'];
resdir = [DATAPATH 'Raphe\zs7\pixelcount_1pc_10px_b\'];
resdir_norm = [DATAPATH 'Raphe\zs7\normalization_1pc_10px_b\'];
mm = pwd;
dr = dir(inpdir_gfp);
dr = dr(3:end);
sf = length(dr);

% Load relative pyramidal layer postition
ff = [inpdir_relpyrpos 'zs7g_' cnm '_RPP.mat'];
load(ff)
mpop = floor(mpostpyr);
mprp = floor(mprepyr);
lm = mpop + mprp + 1;

% Load images
asIg = zeros(sf,lm);
asIr = zeros(sf,lm);
asoIr = zeros(sf,lm);
rtI = zeros(sf,lm);
for o = 1:sf
    fn_gfp = [inpdir_gfp dr(o).name];
    fn_ser = [inpdir_ser dr(o).name(1:end-7) 'ser'];
    [Igreen,cmap] = imread(fn_gfp,'tif');
    [Ired,cmap] = imread(fn_ser,'tif');
    Igreen = double(Igreen);
    Ired = double(Ired);
    Igreen = Igreen / max(Igreen(:));
    Ired = Ired / max(Ired(:));
    
% Load pyramidal layer position
    [pth fnm_gfp xtn] = fileparts(fn_gfp);
    [pth fnm_ser xtn] = fileparts(fn_ser);
    fn_pp = [inpdir_pyrpos fnm_gfp '_PP.mat'];
    load(fn_pp)
    pl = round(pyrpos);

% Remove background and count pixels
    [Ig Igo] = processcolor(Igreen,2);    % remove background
    [Ir Iro] = processcolor(Ired,1);
    fn_n = [resdir_norm fnm_ser '_norm.tif'];
    imwrite(Ir,fn_n,'tif')
    fn_n = [resdir_norm fnm_gfp '_norm.tif'];
    imwrite(Ig,fn_n)
    
    oIgo = Igo;     % overlay
    oIgo(Iro==0) = 0;
    oIro = Iro;
    oIro(Igo==0) = 0;
    Io = zeros(size(Ired));
    Io(:,:,1) = oIro;
    Io(:,:,2) = oIgo;
    
    sIg = sum(Igo);
    sIg_prepl = sIg(pl-mprp:pl);    % below the pyramidal layer
    sIg_postpl = sIg(pl+1:pl+mpop);  % above pyramidal layer
    asIg(o,:) = [sIg_prepl sIg_postpl];
    sIr = sum(Iro);
    sIr_prepl = sIr(pl-mprp:pl);    % below the pyramidal layer
    sIr_postpl = sIr(pl+1:pl+mpop);  % above pyramidal layer
    asIr(o,:) = [sIr_prepl sIr_postpl];
    soIr = sum(oIro);
    soIr_prepl = soIr(pl-mprp:pl);    % below the pyramidal layer
    soIr_postpl = soIr(pl+1:pl+mpop);  % above pyramidal layer
    asoIr(o,:) = [soIr_prepl soIr_postpl];
    rtI(o,:) = asoIr(o,:) ./ asIg(o,:);
end

% Plot result
masIg = mean(asIg);
smasIg = smooth(masIg,'linear');
masIr = mean(asIr);
smasIr = smooth(masIr,'linear')
masoIr = mean(asoIr);
smasoIr = smooth(masoIr,'linear')
mrtI = nanmean(rtI);
smrtI = smooth(mrtI,'linear');
cd(resdir)
dbclear if error
figure
hold on
plot(((1:lm)-mprp)*0.31,smasIg,'g')
plot(((1:lm)-mprp)*0.31,smasIr,'r')
plot(((1:lm)-mprp)*0.31,smasoIr,'y')
saveas(gcf,['zs7g_' cnm '_all.fig'])
saveas(gcf,['zs7g_' cnm '_all.tif'])
figure
plot(((1:lm)-mprp)*0.31,smrtI,'k')
saveas(gcf,['zs7g_' cnm '_ypg.fig'])
saveas(gcf,['zs7g_' cnm '_ypg.tif'])



% ------------------------------------------------------------------------
function [I Ig] = processcolor(I,ch)

% Remove background
Ig = I(:,:,ch);
sig = sort(Ig(:));
lig = length(Ig(:));
levs = [0.95 0.99 0.995 0.999];  % significance levels
lle = length(levs);
sl = zeros(lle,2);
for k = 1:lle
    sl(k,1) = levs(k);
    inx = ceil(lig*levs(k));
    sl(k,2) = sig(inx);     % critical values
end

Ig(Ig<sl(2,2)) = 0;     % thresholding
Ig(Ig>=sl(2,2)) = 1;

Ig = bwareaopen(Ig,11);  % subtraction (remove by area); remove 1 (<2) px elements

I(:,:,ch) = Ig;



% -------------------------------------------------------------------------
function [X2 S] = smooth(X,str)

n = 65;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
S = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
            S(k) = std(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
                S(k) = std([X(mod2(k-nn,m):m); X(1:k+nn)]) / sqrt(n);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
                S(k) = std([X(k-nn:m); X(1:mod2(k+nn,m))]) / sqrt(n);
            else
                X2(k) = mean(X(k-nn:k+nn));
                S(k) = std(X(k-nn:k+nn)) / sqrt(n);
            end
        end
end