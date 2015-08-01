function somresponseprofiles(cellids)
%SOMRESPONSEPROFILES   Overall response profile for various events.
%   This analysis uses an event-spike cross correlation measure to
%   calculate selectivity of neurons to various events in the task.
%
%   Cross-correlations (CCG) are calculated for 5 events:
%   1. HomeZoneOut
%   2. TriggerZoneOut
%   3. HomeZoneIn
%   4. WaterValveOn
%   5. TriggerZoneIn
%   Window size of +-0.5s and bin width of 50 ms is used. The CCGs are
%   plotted and saved for all cells.
%
%   Significance of modulation is assessed by crossing 0.01 significance levels
%   calculated based on shifted crosscorrelations (shifts ranged from 10 to
%   30s). A CCG strength matrix is calculated and saved, containing the following
%   measures (3rd dim.) for all cells (1st dim.) and events (2nd dim.): 
%   1. CCG area outside the significance limits
%   2. CCG area above the upper sign. level
%   3. CCG area below the lower sign. level
%   4. SD of CCG
%   5. sum of CCG
%   The sum of all 5 CCGs for each cell is also calculated and stored.
%
%   See also SOMRESPNSECCG.

% Directories
global DATADIR
if nargin < 1   % load cellIDs
    load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
%     cellids = pv_beh;
%     cellids = som_beh([1:19 21:end]);
    cellids  = listtag('cells');
end
global DATAPATH
resdir = [DATAPATH 'SOM\ResponseCCG_temp\'];    % results directory

% Decide the trigger events, and previous and nextevent - ONLY FIXED EVENTS!
iE = 1;
events{iE} = {'HomeZoneOut1' 'HomeZoneIn1' 'NextTriggerZoneIn'};
iE = iE + 1;
events{iE} = {'TriggerZoneOut' 'TriggerZoneIn' 'HomeZoneIn1'};
iE = iE + 1;
events{iE} = {'HomeZoneIn1' 'TriggerZoneOut' 'HomeZoneOut1'};
iE = iE + 1;
events{iE} = {'WaterValveOn' 'TriggerZoneOut' 'HomeZoneOut1'};
iE = iE + 1;
events{iE} = {'TriggerZoneIn' 'PreviousHomeZoneOut' 'TriggerZoneOut'};

% Arguments
win = [-2 2];
g.eventtype = 'behav';
g.sigma = 0;
g.window = win;
win_margin = [0 0];
g.dt = 0.01;
margin = g.sigma*3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array

% Cell loop
NumCells = length(cellids);
NumEvents = length(events);
CcgStrength = nan(NumCells,NumEvents,4);
SumCcgs = nan(NumCells,1);
problemID = {};
errorstack = {};
for iCell = 1:NumCells
    cellid = cellids(iCell);
    try
        [uapluslas uas las stds sccr sumccrs] = main(cellid,g,events,NumCells,NumEvents,resdir);
        CcgStrength(iCell,:,1) = uapluslas;     % area of CCG outside significance levels
        CcgStrength(iCell,:,2) = uas;   % area of CCG above significance level
        CcgStrength(iCell,:,3) = las;   % area of CCG below significance level
        CcgStrength(iCell,:,4) = stds;     % measure for event-related modulation of CCG
        CcgStrength(iCell,:,5) = sccr;     % sum of individual CCGs
        SumCcgs(iCell) = sumccrs;   % sum of CCGs for future normalization
        
        % Save figure
        cellidt = strrep(cellid{1},'.','_');
        fnm = [resdir cellidt '_RESPONSE_CCG.fig'];
        hgsave2(gcf,fnm)
    catch ME
        disp(ME.message)
        problemID = [problemID cellid];
        errorstack = [errorstack {ME}];
    end
    
    % Save CCG strength measure
    fnm = [resdir 'CCG_STRENGTH.mat'];
    save(fnm,'CcgStrength','SumCcgs')
    fnm = [resdir 'PROBLEMID.mat'];
    save(fnm,'problemID','errorstack')
end

% Plot summary
UaStrength = squeeze(CcgStrength(:,:,4));
mx = max(UaStrength,[],2);
for k = 1:NumCells
    if ~isequal(sum(UaStrength(k,:)),0)
        UaStrength(k,:) = UaStrength(k,:)==mx(k);
    end
end
figure;
plot(nansum(UaStrength))

% -------------------------------------------------------------------------
function [uapluslas uas las stds sccr sumccrs] = main(cellid,g,events,NumCells,NumEvents,resdir)

% Load prealigned spikes
switch g.eventtype
    case 'stim'
        TE=loadcb(cellid,'StimEvents');
        SP=loadcb(cellid,'STIMSPIKES');
    case {'event','behav'}
        TE=loadcb(cellid,'TrialEvents');
        SP=loadcb(cellid,'EVENTSPIKES');
end
figure(39)
clf
setmyfigure(gcf)
sub_h = set_subplots(2,ceil(length(events)/2),0.1,0.1);

% Event loop
uapluslas = nan(1,NumEvents);
uas = nan(1,NumEvents);
las = nan(1,NumEvents);
stds = nan(1,NumEvents);
sccr = nan(1,NumEvents);
sumccrs = 0;
for iE = 1:NumEvents
    
    % Get valid windows for Previous and Next Events
    TriggerEvent = events{iE}{1};
    
    % Get valid trials based on windows
    % wierd isnan accounts for either window being a NaN
    % if Previous event came after triggerevent and
    % if next event came before triggerevent
    valid_i = ~isnan(TE.(TriggerEvent));
%     stimes  = SP.event_stimes{trigger_pos}(valid_i);
%     ev_windows = ev_windows(valid_i,:);
%         binraster = stimes2binraster(stimes,time,g.dt);
    
    % Time series for cross-correlation
    evtimes = TE.(TriggerEvent)(valid_i);
%     stimes2 = arrayfun(@(s)stimes{s}+evtimes(s),1:length(stimes),'UniformOutput',false);
%     stimes2 = cell2mat(stimes2');
    stimes2 = loadcb(cellid,'SPIKES');
    
    % Cross-correlation
    [H1 ccr lwr upr rccg] = somccg_conf_filter(evtimes',stimes2,500,10000,30000);
    
    % Cross-correlation strength measures
    pua1 = ccr - upr;
    pua2 = pua1(pua1>0);
    ua = sum(pua2);   % area of CCG above upper significance level
    pla1 = lwr - ccr;
    pla2 = pla1(pla1>0);
    la = sum(pla2);   % area of CCG below lower significance level
    uapluslas(iE) = ua + la;       % area of CCG outside significance levels
    uas(iE) = ua;
    las(iE) = la;
    stds(iE) = std(ccr);     % measure for event-related modulation of CCG
    sccr(iE) = sum(ccr);     % sum of CCG
    sumccrs = sumccrs + sum(ccr);   % sum of CCGs for future normalization
    
    % Put SD as text on the subplots
    y_lim = ylim;
    text(0,y_lim(2)*0.75,num2str(stds(iE)),'Color','red')
    
    % Plot CCG
    copyobj(allchild(gca),sub_h(iE))
    close(H1)
    axes(sub_h(iE))
    axis tight
    title(events{iE}(1));
    
    % Save CCG
    event = events{iE}{1};
    cellidt = strrep(cellid{1},'.','_');
    fnm = [resdir cellidt '_' event '_CCG.mat'];
    save(fnm,'ccr','lwr','upr','rccg','cellid','event')
end

% -------------------------------------------------------------------------
function [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn,shuff_min,shuff_max,cutoff)
%SOMCCG_CONF_FILTER   Crosscorrelation.
%   SOMCCG_CONF_FILTER(VD1,VD2) calculates crosscorrelogram for
%   discriminated units VD1 and VD2, using a +-50 ms time window.
%   Confidance intervals are calculated based on crosscorrelations of
%   shifted data. SOMCCG_CONF_FILTER calculates confidence interval from
%   high-pass filtered CCG, adding the low-pass filtered part back
%   afterwards (see XCORR_WRAND_FILTER). Filtering is at CUTOFF Hz (see
%   below).
%
%   H = SOMCCG_CONF_FILTER(VD1,VD2) returns the handles of the resulting
%   plot.
%
%   H = SOMCCG_CONF_FILTER(VD1,VD2,WN) uses WN input argument as time
%   window for the cross-correlogram. WN should be given in milliseconds.
%   Normalized crosscorrelogram will not be effected.
%
%   [H1 CCR LWR UPR RCCG] = SOMCCG_CONF_FILTER(VD1,VD2,WN) returns
%   cross-correlogram (CCR), lower and upper confidance intervals (LWR and
%   UPR) and shifted cross-correlations (RCCG).
%
%   [H1 CCR LWR UPR RCCG] = SOMCCG_CONF_FILTER(VD1,VD2,WN,SHUFF_MIN,SHUFF_MAX) 
%   accepts input arguments for minimal and maximal shifts for shuffled
%   crosscorrelations in ms (default: 100 ms and 5 s).
%
%   [H1 CCR LWR UPR RCCG] = SOMCCG_CONF_FILTER(VD1,VD2,WN,SHUFF_MIN,SHUFF_MAX,CUTOFF)
%   accepts a CUTOFF parameter for filtering CCGs (see XCORR_WRAND_FILTER;
%   default is 4 Hz).
%
%   See also SOMCCG, XCORR_WRAND_FILTER and XCORR.

% Input argument check
error(nargchk(2,6,nargin))
if nargin < 3
    wn = 50;    % window size in ms
end
if nargin < 4
    shuff_min = 100;     % 100 ms
end
if nargin < 5
    shuff_max = 5000;    % 5 s
end
if nargin < 6
    cutoff = 4;     % default filter cutoff: 4 Hz
end
sr = 1000;
shuff_min = shuff_min * sr / 1000;
shuff_max = shuff_max * sr / 1000;
num_shuff = 2000;
bns = 50 * sr / 1000;    % bin size for CCG (50 ms)

% Calculate spike times in milliseconds
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;
mn = min(nc1(1),nc2(1));  % only relative spike times count; avoid out of memory
nc1 = nc1 - mn;
nc2 = nc2 - mn;
nc1(nc1<0.5) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
nc2(nc2<0.5) = [];
wn2 = wn / 1000;    % window size in seconds

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
shf = round(rand(1,num_shuff)*shuff_max+shuff_min);
% [ccr,lags,rccg] = xcorr_wrand_filter(sr,cutoff,zunit2,zunit1,wn2*sr,shf);     % 1->2; window: -wn ms - wn ms
[ccr,lags,rccg] = xcorr_wrand(zunit2,zunit1,wn2*sr,shf);
% ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.
ccr = reshape(ccr(1:end-1),bns,(length(ccr)-1)/bns);     % 10 ms bins
ccr = sum(ccr);
rccg = reshape(rccg(:,1:end-1),num_shuff,bns,(size(rccg,2)-1)/bns);
rccg = squeeze(sum(rccg,2));

% Plot
H1 = figure;
time = linspace(-wn,wn,length(ccr));
bar(time,ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])

% Plot confidence interval
hold on
[lc nc] = size(rccg);
ptc = ceil(lc*0.005);
upr = zeros(1,nc);
lwr = zeros(1,nc);
for k = 1:nc
    sts = sort(rccg(:,k),'ascend');
    upr(k) = sts(end-ptc);
    lwr(k) = sts(ptc);
end
plot(time,upr,'Color',[0.7 0.7 0.7])
plot(time,lwr,'Color',[0.7 0.7 0.7])

% -------------------------------------------------------------------------
function [c,lags,rccg] = xcorr_wrand_filter(sr,cutoff,x,varargin)
%XCORR_WRAND_FILTER Cross-correlation function estimates.
%   C = XCORR_WRAND_FILTER(A,B), where A and B are length M vectors (M>1),
%   returns the length 2*M-1 cross-correlation sequence C. If A and B are
%   of different length, the shortest one is zero-padded.  C will be a row
%   vector if A is a row vector, and a column vector if A is a column
%   vector.
%
%   XCORR_WRAND_FILTER produces an estimate of the correlation between two
%   random (jointly stationary) sequences:
%          C(m) = E[A(n+m)*conj(B(n))] = E[A(n)*conj(B(n-m))]
%   It is also the deterministic correlation between two deterministic
%   signals.
%
%   XCORR_WRAND_FILTER(A), when A is a vector, is the auto-correlation
%   sequence. XCORR_WRAND_FILTER(A), when A is an M-by-N matrix, is a large
%   matrix with 2*M-1 rows whose N^2 columns contain the cross-correlation
%   sequences for all combinations of the columns of A. The zeroth lag of
%   the output correlation is in the middle of the sequence, at element or
%   row M.
%
%   XCORR_WRAND_FILTER(...,MAXLAG) computes the (auto/cross) correlation
%   over the range of lags:  -MAXLAG to MAXLAG, i.e., 2*MAXLAG+1 lags. If
%   missing, default is MAXLAG = M-1.
%
%   [C,LAGS] = XCORR_WRAND(...)  returns a vector of lag indices (LAGS).
%
%   XCORR_WRAND_FILTER(...,SCALEOPT), normalizes the correlation according
%   to SCALEOPT:
%     'biased'   - scales the raw cross-correlation by 1/M. 'unbiased' -
%     scales the raw correlation by 1/(M-abs(lags)). 'coeff'    -
%     normalizes the sequence so that the auto-correlations
%                  at zero lag are identically 1.0.
%     'none'     - no scaling (this is the default).
%
%   [C,LAGS,RCCG] = XCORR_WRAND_FILTER(SR,CUTOFF,X,Y,MAXLAG,SCALEOPT,RANDVEC) 
%   or XCORR_WRAND_FILTER(SR,CUTOFF,X,Y,MAXLAG,RANDVEC) uses RANDVEC to
%   calculate shifted cross-correlations (RCCG), which can be used for
%   confidence interval calculations. RANDVEC should contain numbers
%   greater than MAXLAG but less than 10000. The rows of RCCG correspond to
%   the elements of RANDVEC. Shifted cross-correlations are calculated from
%   full cross-correlogram high-pass filtered at CUTOFF [Hz] (to remove
%   drifts and comodulation by slow oscillations; recommendation for
%   short-range CCGs: 4 Hz). The low-pass filtered cross-correlogram is
%   added back to the confidence limits. The sampling rate should be given
%   in SR.
%
%   See also XCORR, XCOV, CORRCOEF, CONV, CCONV, COV and XCORR2.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.16.4.6 $  $Date: 2009/07/14 04:00:10 $ 
%
%   References:
%     S.J. Orfanidis, "Optimum Signal Processing. An Introduction"
%     2nd Ed. Macmillan, 1988.

error(nargchk(1,6,nargin,'struct'));

[x,nshift] = shiftdim(x);
[xIsMatrix,autoFlag,maxlag,scaleType,randVec,msg] = parseinput(x,varargin{:});
if ~isempty(msg), error(generatemsgid('SigErr'),msg); end

if xIsMatrix,
	[c,M,N] = matrixCorr(x);
else
	[c,M,N] = vectorXcorr(x,autoFlag,varargin{:});
end

% Force correlation to be real when inputs are real
c = forceRealCorr(c,x,autoFlag,varargin{:});
c2 = [c(end-10000+1:end,:); c(1:10000+1,:)];
nqf = sr / 2;   % Nyquist frequency
b1 = fir1(2^12,cutoff/nqf,'low');
cbelowX = filtfilt(b1,1,c2);
b2 = fir1(2^12,cutoff/nqf,'high');
caboveX = filtfilt(b2,1,c2);
mpoint = (length(cbelowX) + 1) / 2;

lags = -maxlag:maxlag;
ll = 2 * maxlag + 1;

% Shifted CCG's for confidence interval
ncr = length(randVec);
rccg = nan(ncr,ll);
for k = 1:ncr
    rccg(k,:) = cbelowX(mpoint-maxlag:mpoint+maxlag) + caboveX(randVec(k)-maxlag:randVec(k)+maxlag);
end

% Keep only the lags we want and move negative lags before positive lags
if maxlag >= M,
	c = [zeros(maxlag-M+1,N^2);c(end-M+2:end,:);c(1:M,:);zeros(maxlag-M+1,N^2)];
else
	c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
end

% Scale as specified
[c,msg] = scaleXcorr(c,xIsMatrix,scaleType,autoFlag,M,maxlag,lags,x,varargin{:});
if ~isempty(msg), error(generatemsgid('SigErr'),msg); end

% If first vector is a row, return a row
c = shiftdim(c,-nshift);

%----------------------------------------------------------------
function [c,M,N] = matrixCorr(x)
% Compute all possible auto- and cross-correlations for a matrix input

[M,N] = size(x);

X = fft(x,2^nextpow2(2*M-1));

Xc = conj(X);

[MX,NX] = size(X);
C = zeros(MX,NX*NX);
for n =1:N,
    C(:,(((n-1)*N)+1):(n*N)) = repmat(X(:,n),1,N).*Xc;
end

c = ifft(C);

%----------------------------------------------------------------
function [c,M,N] = vectorXcorr(x,autoFlag,varargin)
% Compute auto- or cross-correlation for vector inputs

x = x(:);

[M,N] = size(x);

if autoFlag,
	% Autocorrelation
	% Compute correlation via FFT
	X = fft(x,2^nextpow2(2*M-1));
	c = ifft(abs(X).^2);
	
else
	% xcorrelation
	y = varargin{1};
	y = y(:);
	L = length(y);
	
	% Cache the length(x)
	Mcached = M;
	
	% Recompute length(x) in case length(y) > length(x)
	M = max(Mcached,L);
	
    if (L ~= Mcached) && any([L./Mcached, Mcached./L] > 10),

        % Vector sizes differ by a factor greater than 10,
		% fftfilt is faster
		neg_c = conj(fftfilt(conj(x),flipud(y))); % negative lags
		pos_c = flipud(fftfilt(conj(y),flipud(x))); % positive lags
		
		% Make them of almost equal length (remove zero-th lag from neg)
		lneg = length(neg_c); lpos = length(pos_c);
		neg_c = [zeros(lpos-lneg,1);neg_c(1:end-1)];
		pos_c = [pos_c;zeros(lneg-lpos,1)];
		
		c = [pos_c;neg_c];
		
	else
		if L ~= Mcached,
			% Force equal lengths
			if L > Mcached
				x = [x;zeros(L-Mcached,1)];
				
			else
				y = [y;zeros(Mcached-L,1)];
			end									
		end
		
		% Transform both vectors
		X = fft(x,2^nextpow2(2*M-1));
		Y = fft(y,2^nextpow2(2*M-1));
		
		% Compute cross-correlation
		c = ifft(X.*conj(Y));
	end
end
	
%----------------------------------------------------------------
function [c,msg] = scaleXcorr(c,xIsMatrix,scaleType,autoFlag,...
	M,maxlag,lags,x,varargin)
% Scale correlation as specified

msg = '';

switch scaleType,
case 'none',
	return
case 'biased', 
	% Scales the raw cross-correlation by 1/M.
	c = c./M;
case 'unbiased', 
	% Scales the raw correlation by 1/(M-abs(lags)).
	scale = (M-abs(lags)).';
	scale(scale<=0)=1; % avoid divide by zero, when correlation is zero
	
	if xIsMatrix,
		scale = repmat(scale,1,size(c,2));
	end
	c = c./scale;
case 'coeff',
	% Normalizes the sequence so that the auto-correlations
	% at zero lag are identically 1.0.
	if ~autoFlag,
		% xcorr(x,y)
		% Compute autocorrelations at zero lag
		cxx0 = sum(abs(x).^2);
		cyy0 = sum(abs(varargin{1}).^2);
		scale = sqrt(cxx0*cyy0);
		c = c./scale;
	else
		if ~xIsMatrix,
			% Autocorrelation case, simply normalize by c[0]
			c = c./c(maxlag+1);
		else
			% Compute the indices corresponding to the columns for which
			% we have autocorrelations (e.g. if c = n by 9, the autocorrelations
			% are at columns [1,5,9] the other columns are cross-correlations).
			[mc,nc] = size(c);
			jkl = reshape(1:nc,sqrt(nc),sqrt(nc))';
			tmp = sqrt(c(maxlag+1,diag(jkl)));
			tmp = tmp(:)*tmp; 
			cdiv = repmat(tmp(:).',mc,1);
			c = c ./ cdiv; % The autocorrelations at zero-lag are normalized to
			% one
		end
	end
end	

%----------------------------------------------------------------
function [xIsMatrix,autoFlag,maxlag,scaleType,randVec,msg] = parseinput(x,varargin)
%    Parse the input and determine optional parameters:
%
%    Outputs:
%       xIsMatrix - flag indicating if x is a matrix
%       AUTOFLAG  - 1 if autocorrelation, 0 if xcorrelation
%       maxlag    - Number or lags to compute
%       scaleType - String with the type of scaling wanted
%       randVec   - array of random shifts for conf. interval
%       msg       - possible error message

% Set some defaults
scaleType = '';
autoFlag = 1; % Assume autocorrelation until proven otherwise
maxlag = []; 
randVec = [];
xIsMatrix = false;

errMsg = 'Input argument is not recognized.';

switch nargin,
case 2,
	% Can be (x,y), (x,maxlag), or (x,scaleType)
	if ischar(varargin{1}),
		% Second arg is scaleType
		scaleType = varargin{1};
		
	elseif isnumeric(varargin{1}),
		% Can be y or maxlag
		if length(varargin{1}) == 1,
			maxlag = varargin{1};
		else
			autoFlag = 0;
			y = varargin{1};
		end
	else
		% Not recognized
		msg = errMsg;
		return
	end
case 3,
	% Can be (x,y,maxlag), (x,maxlag,scaleType) or (x,y,scaleType)
	maxlagflag = 0; % By default, assume 3rd arg is not maxlag
	if ischar(varargin{2}),
		% Must be scaletype
		scaleType = varargin{2};
		
	elseif isnumeric(varargin{2}),
		% Must be maxlag
		maxlagflag = 1;
		maxlag = varargin{2};
		
	else
		% Not recognized
		msg = errMsg;
		return
	end
	
	if isnumeric(varargin{1}),
		if maxlagflag,
			autoFlag = 0;
			y = varargin{1};
		else
			% Can be y or maxlag
			if length(varargin{1}) == 1,
				maxlag = varargin{1};
			else
				autoFlag = 0;
				y = varargin{1};
			end
		end
	else
		% Not recognized
		msg = errMsg;
		return
	end
	
case 4,
	% Can be (x,y,maxlag,scaleType) or (x,y,maxlag,randVec)
	autoFlag = 0;
	y = varargin{1};
	
	maxlag = varargin{2};
	
    if ~isnumeric(varargin{3})
        scaleType = varargin{3};
    else
        randVec = varargin{3};
    end
    
case 5
    % Must be (x,y,maxlag,scaleType,randVec)
    autoFlag = 0;
	y = varargin{1};
	
	maxlag = varargin{2};
	
    scaleType = varargin{3};
    
    randVec = varargin{4};
end

% Determine if x is a matrix or a vector
[xIsMatrix,m] = parse_x(x);



if ~autoFlag,
	% Test y for correctness
	[maxlag,msg] = parse_y(y,m,xIsMatrix,maxlag);
	if ~isempty(msg),
		return
	end
end

[maxlag,msg] = parse_maxlag(maxlag,m);
if ~isempty(msg),
	return
end


% Test the scaleType validity
[scaleType,msg] = parse_scaleType(scaleType,errMsg,autoFlag,m,varargin{:});
if ~isempty(msg),
	return
end

	
%-------------------------------------------------------------------
function [xIsMatrix,m] = parse_x(x)


xIsMatrix = (size(x,2) > 1);

m = size(x,1);


%-------------------------------------------------------------------
function [maxlag,msg] = parse_y(y,m,xIsMatrix,maxlag)
msg = '';
[my,ny] = size(y);
if ~any([my,ny] == 1),
	% Second arg is a matrix, error
	msg = 'B must be a vector (min(size(B))==1).';
	return
end

if xIsMatrix,
	% Can't do xcorr(matrix,vector)
	msg = 'When B is a vector, A must be a vector.';
	return
end
if (length(y) > m) && isempty(maxlag),
	% Compute the default maxlag based on the length of y
	maxlag = length(y) - 1;
end

%-------------------------------------------------------------------
function [maxlag,msg] = parse_maxlag(maxlag,m)
msg = '';
if isempty(maxlag),
	% Default hasn't been assigned yet, do so
	maxlag = m-1;
else	
	% Test maxlag for correctness	
	if  length(maxlag)>1
		msg = 'Maximum lag must be a scalar.';
		return
	end	
	if maxlag < 0,
		maxlag = abs(maxlag);
	end
	if maxlag ~= round(maxlag),
		msg = 'Maximum lag must be an integer.';
	end
end

%--------------------------------------------------------------------
function c = forceRealCorr(c,x,autoFlag,varargin)
% Force correlation to be real when inputs are real

forceReal = 0; % Flag to determine whether we should force real corr

if (isreal(x) && autoFlag) || (isreal(x) && isreal(varargin{1})),
	forceReal = 1;
end


if forceReal,
	c = real(c);
end

%--------------------------------------------------------------------
function [scaleType,msg] = parse_scaleType(scaleType,errMsg,autoFlag,m,varargin)
msg = '';
if isempty(scaleType),
	scaleType = 'none';
else
	scaleOpts = {'biased','unbiased','coeff','none'};
	indx = find(strncmpi(scaleType, scaleOpts, length(scaleType)));
	
	if isempty(indx),
		msg = errMsg;
		return
	else
		scaleType = scaleOpts{indx};
	end
	
	if ~autoFlag && ~strcmpi(scaleType,'none') && (m ~= length(varargin{1})),
		msg = 'Scale option must be ''none'' for different length vectors A and B.';
		return
	end
end