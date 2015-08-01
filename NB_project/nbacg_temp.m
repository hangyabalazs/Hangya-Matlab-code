function nbacg(cellids,varargin)
%NBACG   Auto-correlation.
%   NBACG calculates cross-correlations. For details on the algorithm, 
%   see XCORR.
%
%   NBACG(I) calls ACG for cells defined by the cell ID list (or index set
%   to CELLIDLIST, see CellBase documentation) I.
%   Optional input arguments (parameter-value pairs with default values):
%       'issave', false - controls saving behavior, plots and
%           auto-correlation matrices are saved only if 'issave' is set to
%           true
%
%   See also NBCCG and XCORR.

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addParamValue(prs,'issave',false,@islogical)   % control saving behavior
parse(prs,cellids,varargin{:})
g = prs.Results;

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB\ACG\medium\'];  % results directory
fnmm = 'ACG_matrices.mat';   % filename for saving the result matrices

% Include only long enough segments
longsegments = false;  % control whether to use this option
seglim = 300;

% Determine time window
sr = 1000;      % sampling rate
wn = 500 * sr / 1000;    % 2 * 50 ms window
res = 0.5;   % resolution for ACG in ms

% Input argument check
if nargin < 1
    loadcb   % load CellBase
    cellids = CELLIDLIST; 
else
    if isnumeric(cellids)
        loadcb   % load CellBase
        cellids = CELLIDLIST(cellids);
    end
end
if ischar(cellids)
    cellids = {cellids};  % one cell ID
end

% Cell loop for ACG
wb = waitbar(0,'Please wait...','Name','Running NBACG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [1 50000];   % include max 50000 spikes; no lower limit
numCells = length(cellids);
[CCR SCCR] = deal(zeros(numCells,2*wn/res));
[SegmentLength BurstIndex Refractory ThetaIndex] = deal(nan(numCells,1));
for iC = 1:numCells   % loop through the cells
    cell = cellids{iC};
    try
        tseg = findSegs3(cell,'segfilter','stim_excl_nb',...
            'light_activation_duration',[-5 5],'margins',[0 0]);  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if longsegments
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < seglim   % use the longest segment if it's longer than the threshold
                continue
            end
        end
        SegmentLength(iC) = sum(ltseg);  % cumulative length of the segments
        [ncc,seltsind,selisi] = extractSegSpikes(cell,tseg);   % find spikes in the time segments
    catch ME
        disp(ME.message)
        disp('Could not extract the segment.');
    end

%     ncc = loadcb(cell,'SPIKES');   % use all spikes
        
    if length(ncc) > limit_spikes(2);      % crop if too long to avoid out of memory
        ncc = ncc(1:limit_spikes(2));
    end
    
    if length(ncc) > limit_spikes(1)     % minimum criterion
        [H1 ccr lags] = acorr(ncc,wn,res);
        sccr = smooth(ccr,'linear',21);    % smoothed ACR
%         nqf = 1 / res * 1000 / 2;   % Nyquist freq.
%         flt = fir1(32,[4 10]/nqf,'bandpass');
%         sccr2 = filtfilt(flt,1,sccr);   % high-pass filter > 1 Hz
        hold on
        plot(lags,sccr,'Color',[0.7 0.7 0.7])
%         plot(lags,sccr2,'c')
        bar(lags(lags>-10&lags<10),ccr(lags>-10&lags<10),'FaceColor','g','EdgeColor','g')
        bar(lags(lags>-200&lags<-180),ccr(lags>-200&lags<-180),'FaceColor','r','EdgeColor','r')
        bar(lags(lags>180&lags<200),ccr(lags>180&lags<200),'FaceColor','r','EdgeColor','r')
        bar(lags((lags>=180&lags<=200)|(lags>=65&lags<=85)),...
            ccr((lags>=180&lags<=200)|(lags>=65&lags<=85)),...
            'FaceColor','c','EdgeColor','c');
        BurstIndex(iC) = burstinx(ccr,lags);   % burst index
        Refractory(iC) = refract(ccr,sccr,lags);   % refractory
        ThetaIndex(iC) = thetainx4(sccr,lags);   % theta index
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.95,regexprep(cell,'_',' '))
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,['Burst index: ' num2str(BurstIndex(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,['Refractory: ' num2str(Refractory(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.8,['Theta index: ' num2str(ThetaIndex(iC))])
        if g.issave   % save figure
            ncl = regexprep(cell,'\.','_');
            fnm = ['ACG_' ncl '.fig'];
            saveas(H1,fullfile(resdir,fnm))   % save CCG plot
            close(H1)
        end
        CCR(iC,:) = ccr;   % cross-correlogram
        SCCR(iC,:) = sccr;   % cross-correlogram
        disp(['Cell #' num2str(iC) ' / ' num2str(numCells) ' done......'])
    end
    
    % Save
    if g.issave
        if isequal(mod(iC,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(resdir,fnm),'cellids','CCR','SCCR','lags',...
                'SegmentLength','BurstIndex','Refractory','ThetaIndex')
            disp('Autosave done.')
        end
    end
    
    waitbar(iC/numCells)
end
close(wb)   % eliminate progress indicator

% Save
if g.issave
    save(fullfile(resdir,fnmm),'cellids','CCR','SCCR','lags',...
        'SegmentLength','BurstIndex','Refractory','ThetaIndex')
end

% -------------------------------------------------------------------------
function [H1 ccr lags] = acorr(ncc,wn,res)

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;
mn = nc(1);  % only relative spike times count; avoid out of memory
nc = nc - mn;
nc(nc<0.5) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
wn2 = wn / 1000;    % window size in seconds

% Auto-correlogram
zunit = zeros(1,round(nc(end)/res)+5);
zunit(round(nc/res)) = 1;
[ccr lags] = xcorr(zunit,zunit,wn2*sr/res);     % 1->2; window: -wn ms - wn ms
ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
lags(length(lags)/2+0.5) = [];
lags = lags * res;   % in ms

% Plot
H1 = figure;
bar(lags,ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])

% -------------------------------------------------------------------------
function bi = burstinx(ccr,lags)

a = max(ccr(lags>0&lags<10));
b = mean(ccr(lags>180&lags<200));
bi = (a - b) / max(a,b);   % burst index

% -------------------------------------------------------------------------
function tau = refract(ccr,sccr,lags)

mx = max(sccr);   % peak: from smoothed ACG
halfmx = mx / 2;   % half-max
ccr2 = ccr(lags>0);   % half ACG
lags2 = lags(lags>0);   % corresponding lag values
tau = lags2(find(ccr2>halfmx,1,'first'));   % refractory

% -------------------------------------------------------------------------
function thinx = thetainx(ccr,lags)

% Sampling rate
res = lags(2) - lags(1);  % ACG resolution
sr = 1 / res * 1000;   % ACG sampling rate

% Auto-spectrum
% ccr2 = ccr(lags>0);   % half ACG
[y f] = b_fft2(ccr,sr);
% figure
% plot(f(f>=2&f<=15),y(f>=2&f<=15))   % plot auto-spectrum

% Theta index
thinx = sum(y(f>=4&f<=10)) / sum(y(f>=2&f<=50));

% -------------------------------------------------------------------------
function thinx = thetainx2(ccr,lags)

thpeak = ccr(lags>=100&lags<=250);   % ACG first theta peak
thinx = max(thpeak) / min(thpeak);   % theta index

% -------------------------------------------------------------------------
function thinx = thetainx3(ccr,lags,ncc)

% Normalize
% ccr = ccr / length(ncc);

% Theta index
thpeak = ccr(lags>=100&lags<=250);   % ACG first theta peak
[jnk pl] = max(thpeak);
mthp = mean(thpeak);
mthp = max(mthp,1);  % numbers <1 rounded up
ploc = lags(find(lags>=100,1,'first')+pl-1);   % peak location
thtrough = ccr(lags>=ploc&lags<=2*ploc);   % ACG first theta trough
mtht = max(mean(thtrough),1);  % numbers <1 rounded up
thinx = (mthp - mtht) / max(mthp,mtht);   % theta index

% -------------------------------------------------------------------------
function thinx = thetainx4(ccr,lags)

% Theta index
thpeak = ccr(lags>=100&lags<=200);   % ACG first theta peak
[jnk pl] = max(thpeak);
ploc = lags(find(lags>=100,1,'first')+pl-1);   % peak location
mthp = mean(ccr(lags>=ploc-25&lags<ploc+25));

thtrough = ccr((lags>=180&lags<=200)|(lags>=65&lags<=85));   % ACG first theta trough
mtht = mean(thtrough);
thinx = (mthp - mtht) / max(mthp,mtht);   % theta index