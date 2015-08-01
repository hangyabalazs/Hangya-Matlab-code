function bimodalITI6
%BIMODALITI6   Reaction time analysis.
%   BIMODALITI6 analyzes reaction (RT) time in the function of the
%   foreperiod. All sessions with bimodal foreperiod distributions are
%   included; trials are pooled from tese sessions. BIMODALITI6 plots
%   median RT against foreperiod length for a specific stimulus intensity,
%   smoothed with a moving average. Errorshade shows SE of the median,
%   calculated with bootstrap method (1000 resampling).
%
%   See also BIMODALITI7.

% Specify folder
basedir = getpref('cellbase','datapath');
mousename = 'n023';

% Sessions
% Animal 'n023'
SESSIONS.n023.bimodal0 = {'111220a' '111221a' '111222a' '111223a' '111227a' '111228a' ...
    };    % bimodal, different mixing proportions (narrow peaks)
SESSIONS.n023.bimodal = {'111229a' '111230a' '111231a' '120101a' '120102a' '120103a' ...
    '120104a' '120105a' '120106a' '120110a' '120111a' '120112a' '120113a' ...
    '120114a' '120115a' '120116a'};    % bimodal
SESSIONS.n023.unimodal0 = {'120117a' '120118a' '120119a' '120120a' '120121a' '120122a'};   % unimodal
SESSIONS.n023.unimodal = {'120123a' '120124a' '120125a' '120126a' '120127a' '120128a'...
    '120131a' '120201a' '120202a' '120203a'};   % unimodal with uniform
SESSIONS.n023.exponential = {'120207a' '120208a' '120209a' '120210a' '120211a' ...
    '120213a' '120214a' '120215a' '120216a' '120217a'};     % exponential
SESSIONS.n023.unimodalnew = {'120220a' '120221a'};       % unimodal again

% Animal 'n026'
SESSIONS.n026.exponential = {'120207a' '120208a' '120209a' '120210a' '120211b' ...
    '120213a' '120214a' '120215a' '120216a'};     % exponential
SESSIONS.n026.bimodal = {'120218a' '120220a' '120221a' '120222a' '120223a'};  % bimodal (not performing after Cosyne, not included)
        
% Animal 'n027'
SESSIONS.n027.exponential = {'120208a' '120209a' '120210a' '120211a'};     % exponential
SESSIONS.n027.bimodal = {'120213b' '120214a' '120215a' '120215c' '120216a' '120217a' ...
    '120218a' '120220a' '120221a' '120222a' '120223a' ...
    '120301a' '120302a' '120303a' '120305a' '120306a'...
    '120307a' '120308a' '120309a' '120310a' '120311a'...
    '120313a' '120314a' '120315a'};     % bimodal (120215a and c: not doing)
        
% Animal 'n028'
SESSIONS.n028.exponential = {'120207a' '120208a' '120210a' '120211a' ...
    '120213a' '120214a' '120215a' '120216a'};     % exponential
SESSIONS.n028.bimodal = {'120217a' '120218a' '120220a' '120221a'...
    '120222a' '120222c' '120223a' ...
    '120301a' '120302a' '120303a' '120305a' '120306a' '120307a'};     % bimodal
        
% Animal 'n029'
SESSIONS.n029.bimodal = {'120204a' '120204b' '120206a' '120207a' '120207b' ...
    '120208a' '120209a' '120210a' '120211a' '120213a'...
    '120213b' '120214a' '120214b' '120215a' '120216b' '120217a' ...
    '120218a' '120220a' '120220b' '120221a' '120221b' '120222a' '120222b' ...
    '120301a' '120302a' '120303a' '120305a'...
    '120306a' '120307a' '120308a' '120309a' '120310a' '120311a'...
    '120313a' '120314a' '120315a'};     %#ok<STRNU> % bimodal (remark: 120209a session very long and good RT data!!!)

eval(['sessions = SESSIONS.' mousename '.bimodal;']);
numSessions = length(sessions); %#ok<USENS>

% Preallocate
meanrt = cell(numSessions,1);
meanhit = cell(numSessions,1);
meanmiss = cell(numSessions,1);
medianrt = cell(numSessions,1);
medianhit = cell(numSessions,1);
medianmiss = cell(numSessions,1);
nrt = cell(numSessions,1);
nhit = cell(numSessions,1);
nhits = cell(numSessions,1);
nmiss = cell(numSessions,1);
nmisss = cell(numSessions,1);
allrt = cell(numSessions,1);
alliti = cell(numSessions,1);
allhit = cell(numSessions,1);
allmiss = cell(numSessions,1);
allsi = cell(numSessions,1);

allrt2 = [];
alliti2 = [];
allhit2 = [];
allmiss2 = [];
allsi2 = [];

% Binning
edges = 0.1:0.15:5;
cnts = (edges(1:end-1) + edges(2:end)) / 2;
bno = length(cnts);

% Get session data
ints = [20 30 40 50];
clr = summer(length(ints));
clr(2,:) = [0 0 0];
I = 2;  % 30 dB
for iS = 1:numSessions
    session = sessions{iS};
    datapath = fullfile(basedir,mousename,session);
    list = dir(datapath);
    list = setdiff({list.name},{'.';'..';'.DS_Store'});
    stc = strfind(list,'TE.mat');
    cf = ~cellfun(@isempty,stc);
    list = list(cf);
    
    % Load behavioral data
    TE = load(fullfile(datapath,list{1}));
    
    % Find unique sound intensities
    SI = unique(TE.StimulusDuration);
    SIn = SI(~isnan(SI));
    if ~isequal(SIn,ints)
        warning('bimodalITI:differentIntensities','Intensities used differ from default.')
    end
    
    % RT vs ITI
    TE2 = filterTE(TE,TE.StimulusDuration==ints(I));
    TE2 = TE;
    [meanrt{iS} medianrt{iS} nrt{iS} allrt{iS} alliti{iS}] = condmean([TE2.ITIDistribution' TE2.GoRT'],edges);
    [meanhit{iS} medianhit{iS} nhit{iS} allhit{iS} alliti{iS}] = condmean([TE2.ITIDistribution' TE2.Hit'],edges);
    nhits{iS} = cellfun(@(g)sum(nan2zero(g)),allhit{iS});
    [meanmiss{iS} medianmiss{iS} nmiss{iS} allmiss{iS} alliti{iS}] = condmean([TE2.ITIDistribution' TE2.Miss'],edges);
    nmisss{iS} = cellfun(@(g)sum(nan2zero(g)),allmiss{iS});
    [~, ~, ~, allsi{iS}] = condmean([TE2.ITIDistribution' TE2.StimulusDuration'],edges);
    allrt2 = [allrt2 TE2.GoRT];
    alliti2 = [alliti2 TE2.ITIDistribution];
    allhit2 = [allhit2 TE2.Hit];
    allmiss2 = [allmiss2 TE2.Miss];
    allsi2 = [allsi2 TE2.StimulusDuration];
end

% Pool all data bin-wise and calculate median
mnsall = cell(1,bno);
mnsmedian = nan(1,bno);
mnsmedian_se = nan(1,bno);
for t = 1:bno
    for k = 1:numSessions
        mnsall{t} = [mnsall{t}; allrt{k}{t}];
    end
    mnsmedian(t) = nanmedian(mnsall{t});
    if ~isempty(mnsall{t}) && ~all(isnan(mnsall{t}))
        mnsmedian_se(t) = se_of_median(mnsall{t},1000);
    end
end
figure
plot(cnts,mnsmedian)
hold on
errorshade(cnts,mnsmedian,mnsmedian_se,...
    'LineColor',clr(I,:),'ShadeColor',clr(I,:))
[nm xout] = hist(TE.ITIDistribution,50);
stairs(xout,nm/300+0.2,'Color','r')

% Plot smoothed version
smkernel = 3;
smnsmedian = smooth(mnsmedian,'linear',smkernel);
smnsmedian_se = smooth(mnsmedian_se,'linear',smkernel);
figure;
plot(cnts,smnsmedian)
hold on
errorshade(cnts,smnsmedian,smnsmedian_se,...
    'LineColor',clr(I,:),'ShadeColor',clr(I,:))
[nm xout] = hist(TE.ITIDistribution,50);
stairs(xout,nm/900+0.18,'Color','r')

% -------------------------------------------------------------------------
function TE = filterTE(TE,valid_trials)

fnm = fieldnames(TE);
for k = 1:length(fnm)
    TE.(fnm{k}) = TE.(fnm{k})(valid_trials);
end