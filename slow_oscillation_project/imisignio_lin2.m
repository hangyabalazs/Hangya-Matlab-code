function imisignio_lin2(pat,patno,eg,nm_rows,nm_cols,esg)
%IMISIGNIO_LIN2   Input-output function.
%   IMISIGNIO_LIN2 calculates input-output function for significant
%   propagation events based on cross-correlations. Data is plotted and
%   saved for each channels. Background color for multiplot plots is set
%   according to the strength of significant linear correlations. Optional
%   input
%   arguments:
%       PAT: patient code and name
%       PATNO: patient code
%       EG: EEG file no.
%       NM_ROWS: number of grid rows
%       NM_COLS: number of grid columns
%       ESG: segment boundaries
%   Use ICALLER2 for sequential calling!
%
%   See also IMISHIFTCALL_IO and IMISHIFT_IO.

% Input argument check
error(nargchk(0,6,nargin))
if nargin <6
    esg = [1040 1100];
end
if nargin < 5
    nm_cols = 5;           % number of columns on the electrode grid
end
if nargin < 4
    nm_rows = 4;           % number of rows on the electrode grid
end
if nargin < 3
    eg = '46';
end
if nargin < 2
    patno = num2str('40');
end
if nargin < 1
    pat = 'oiti40_mojzsis';
end

% Directories
global DATADIR
global DATAPATH
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\CCGmap\'];
inpdir_lin = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\CCGmap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\IOsign_lin2\'];
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end

% Load MI map
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)

% Significance level
fn = [inpdir 'siglev_EEG' eg];
load(fn)
siglev = sl(4,2);   % sig. lev.: 0.0001

% Transform input variables
Adj = rIMax;        % adjacency matrix
Adj(rIMax<siglev) = NaN;
Adj(isnan(Adj)) = 0;
Adj(Adj>0) = 1;     % time-varying adjacency matrix
rIM = rIMax;
rIM(rIMax<siglev) = NaN;
rIM(isnan(rIM)) = 0;
rIML = rIMaxLoc;
rIML(rIMax<siglev) = NaN;
rIML(isnan(rIML)) = 0;

% Recording properties
chnum = size(Adj,1);   % number of grid electrodes
tmax = size(Adj,3);    % time axis length
if ~isequal(nm_rows*nm_cols,chnum)
    error('Improper number for grid electrodes.')
end
sr = 1000;             % sampling rate

% Load raw data
ddir = [DATADIR '\human_SO\' pat '\grid\mat_EEG_' eg '\'];
fnm = ['EEG_' eg '_' num2str(esg(1)) '_' num2str(esg(2)) '_rs_filt01_40.mat'];
load([ddir fnm])

% Input-output function
cumdatax = cell(chnum,chnum);
cumdatay = cell(chnum,chnum);
for i = 1:chnum
    for j = 1:chnum
        for k = 1:tmax
            if rIM(i,j,k)
                data1 = data((k-1)*100+1:(k-1)*100+1000,i);
                data2 = data((k-1)*100+1+rIML(i,j,k):(k-1)*100+1000+rIML(i,j,k),j);
                cumdatax{i,j} = [cumdatax{i,j} data1];    % for I/O function
                cumdatay{i,j} = [cumdatay{i,j} data2];
            end
        end
    end
end

% Load CCG map
fn = [inpdir_lin 'MIshiftfine.mat'];
load(fn)
fn = [inpdir_lin 'siglev_EEG' eg];
load(fn)
siglev_lin = sl(4,2);   % sig. lev.: 0.0001
rIMlin = rIMax;
rIMlin(rIMax<siglev_lin) = NaN;
rIMlin(isnan(rIMlin)) = 0;

% Plotting and saving I/O function
srIM = squeeze(sum(rIMlin,3));
srIMn = (srIM - min(srIM(:))) / (max(srIM(:)) - min(srIM(:)));
C = colormap('jet');
whitebg(C(10,:))
Hsd = figure;
Hse = figure;
for i = 1:chnum
    for j = 1:chnum
        x = cumdatax{i,j};
        y = cumdatay{i,j};
        mxx = max(x);
        mnx = min(x);
        MC = mnx:((mxx-mnx)/50):mxx;
        next = 1;
        fcn = [];
        fsem = [];
        fsd = [];
        nm = [];
        for mc = MC(2:end);
            xinx = find(x>MC(next)&x<MC(next+1));
            yc = y(xinx);
            fcn(next) = mean(yc);
            fsem(next) = std(yc) / sqrt(length(xinx));
            fsd(next) = std(yc);
            nm(next) = length(xinx);
            next = next + 1;
        end
        fy = (MC(2:end) + MC(1:end-1)) / 2;
        
        spw = 1 / chnum;    % subplot coordinates
        lft = (j - 1) * spw;
        btm = (chnum - i) * spw;
        wth = spw;
        hgt = spw;
        
%         whitebg(C(10,:))
        figure(Hsd)
        set(gcf,'Unit','normalized')
        subplot('Position',[lft btm wth hgt]);
        plot(fy,fcn+fsd,'Color',[0.7 0.7 0.7])
        hold on
        plot(fy,fcn-fsd,'Color',[0.7 0.7 0.7])
        plot(fy,fcn)
        clr = C(min(round(srIMn(i,j)*64)+1,64),:);
        set(gca,'Color',clr)    % set background color to reflect connection strength
        axis tight
        set(gca,'XTick',[],'YTick',[])
%         whitebg(C(10,:))
        figure(Hse)
        set(gcf,'Unit','normalized')
        subplot('Position',[lft btm wth hgt]);
        plot(fy,fcn+fsem,'Color',[0.7 0.7 0.7])
        hold on
        plot(fy,fcn-fsem,'Color',[0.7 0.7 0.7])
        plot(fy,fcn)
        clr = C(min(round(srIMn(i,j)*64)+1,64),:);
        set(gca,'Color',clr)    % set background color to reflect connection strength
        axis tight
        set(gca,'XTick',[],'YTick',[])
    end
end
fnm = 'IOsd.fig';    % save
saveas(Hse,fnm)
fnm = 'IOse.fig';
saveas(Hse,fnm)
fnm = 'IO.mat';
save(fnm,'cumdatax','cumdatay')

for i = 1:chnum
    for j = 1:chnum
        x = cumdatax{i,j};
        y = cumdatay{i,j};
        mxx = max(x);
        mnx = min(x);
        MC = mnx:((mxx-mnx)/50):mxx;
        next = 1;
        fcn = [];
        fsem = [];
        fsd = [];
        nm = [];
        for mc = MC(2:end);
            xinx = find(x>MC(next)&x<MC(next+1));
            yc = y(xinx);
            fcn(next) = mean(yc);
            fsem(next) = std(yc) / sqrt(length(xinx));
            fsd(next) = std(yc);
            nm(next) = length(xinx);
            next = next + 1;
        end
        fy = (MC(2:end) + MC(1:end-1)) / 2;
        
        H1 = figure;    % save
        whitebg([0.7 0.7 0.7])
        plot(fy,fcn)
        set(gca,'Color','white')
        hold on
        errorbar(fy,fcn,fsd,'k+')
        H2 = figure;
        plot(fy,fcn)
        set(gca,'Color','white')
        hold on
        errorbar(fy,fcn,fsem,'k+')
        fnm = ['IOsd_' num2str(i) '_' num2str(j) '.fig'];
        saveas(H1,fnm)
        fnm = ['IOse_' num2str(i) '_' num2str(j) '.fig'];
        saveas(H1,fnm)
        fnm = ['IO_' num2str(i) '_' num2str(j) '.mat'];
        save(fnm,'fcn','fsd','fsem','fy')
        close(H1,H2)
    end
end