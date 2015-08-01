function bimodalITI5
% all bimodal sessions included

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
    '120313a' '120314a' '120315a'};     % bimodal (remark: 120209a session very long and good RT data!!!)

sessions = SESSIONS.n023.bimodal;
numSessions = length(sessions);

% Get session data
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

edges = [0.1:0.05:2.7; 0.6:0.05:3.2];     % phase histogram bin limits
cnts = (edges(1,:) + edges(2,:)) / 2;     % phase histogram bin centers
% cnts = 0:0.01:5;
% edges = edge_generator(cnts);
% edges = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95;...
%          0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05];
% cnts = (edges(1,:) + edges(2,:)) / 2;
% edges = [0.1:0.03:0.6 0.6:0.1:1.6;...
%          0.15:0.03:0.65 0.9:0.1:1.9];
% cnts = (edges(1,:) + edges(2,:)) / 2;

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
    
    % RT vs ITI
    maxint = SIn(end);
    if length(SIn) > 1
        maxint2 = SIn(end-1);
    end
    if length(SIn) > 2
        maxint3 = SIn(end-2);
    end
    if length(SIn) > 3
        maxint4 = SIn(end-3);
    end
    TE = filterTE(TE,TE.StimulusDuration==maxint3);
    [meanrt{iS} medianrt{iS} nrt{iS} allrt{iS} alliti{iS}] = condmean([TE.ITIDistribution' TE.GoRT'],edges);
    [meanhit{iS} medianhit{iS} nhit{iS} allhit{iS} alliti{iS}] = condmean([TE.ITIDistribution' TE.Hit'],edges);
    nhits{iS} = cellfun(@(g)sum(nan2zero(g)),allhit{iS});
    [meanmiss{iS} medianmiss{iS} nmiss{iS} allmiss{iS} alliti{iS}] = condmean([TE.ITIDistribution' TE.Miss'],edges);
    nmisss{iS} = cellfun(@(g)sum(nan2zero(g)),allmiss{iS});
    [~, ~, ~, allsi{iS}] = condmean([TE.ITIDistribution' TE.StimulusDuration'],edges);
    allrt2 = [allrt2 TE.GoRT];
    alliti2 = [alliti2 TE.ITIDistribution];
    allhit2 = [allhit2 TE.Hit];
    allmiss2 = [allmiss2 TE.Miss];
    allsi2 = [allsi2 TE.StimulusDuration];
end

mnsall = cell(1,length(cnts));
mnsmedian = nan(1,length(cnts));
for t = 1:length(cnts)
    for k = 1:numSessions
        mnsall{t} = [mnsall{t}; allrt{k}{t}];
%         mnsall{t}(mnsall{t}<0.001) = NaN;
%         mnsall{t}(mnsall{t}>0.7) = NaN;
    end
    mnsmedian(t) = nanmedian(mnsall{t});
end
figure
plot(cnts,mnsmedian)
hold on
[nm xout] = hist(TE.ITIDistribution,50);
stairs(xout,nm/300+0.2,'Color','r')

[nm xnum] = running_average([alliti2' allrt2']);
figure
plot(xnum,nm)

keyboard

% -------------------------------------------------------------------------
function [nm xnum] = running_average(spvr)

bno = 51;
spvr = sortrows(spvr,1);
naninx = isnan(spvr(:,1)) | isnan(spvr(:,2));
spvr(naninx,:) = [];
% mno = floor(size(spvr,1)/bno);
mno = size(spvr,1) - bno;
nm = nan(1,mno);
xnum = nan(1,mno);
for k = 1:mno
    inx1 = k;
    inx2 = inx1 + bno - 1;
    nm(k) = mean(spvr(inx1:inx2,2));
    xnum(k) = mean(spvr(inx1:inx2,1));
end

% -------------------------------------------------------------------------
function [meanrt medianrt nrt allrt alliti] = dhist(spvr,edges)

n = length(edges);
meanrt = zeros(1,n-1);
medianrt = zeros(1,n-1);
nrt = zeros(1,n-1);
allrt = cell(1,n-1);
for k = 2:n
    inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<edges(k));
    meanrt(k-1) = nanmean(spvr(inx,2));
    medianrt(k-1) = nanmedian(spvr(inx,2));
    nrt(k-1) = length(spvr(inx,2));
    allrt{k-1} = spvr(inx,2);
    alliti{k-1} = spvr(inx,1);
end

% -------------------------------------------------------------------------
function TE = filterTE(TE,valid_trials)

fnm = fieldnames(TE);
for k = 1:length(fnm)
    TE.(fnm{k}) = TE.(fnm{k})(valid_trials);
end

% -------------------------------------------------------------------------
function edges = edge_generator(cnts)

% My bimodal failure distribution

NumTrials = 10000;
ITIMin = 0.1;
ITIMax = 3;
mng1 = 0.3;   % parameters for the Gaussians
mng2 = 2;
sdg = 0.15;
pmx1 = 0.35;   % mixing probabilities
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2;

ITIs1 = random('Normal',mng1,sdg,1,NumTrials);
while any(ITIs1>ITIMax) | any(ITIs1<ITIMin)
    inx = ITIs1 > ITIMax  | ITIs1 < ITIMin;
    ITIs1(inx) = random('Normal',mng1,sdg,1,sum(inx));
end
ITIs2 = random('Normal',mng2,sdg,1,NumTrials);
while any(ITIs2>ITIMax) | any(ITIs2<ITIMin)
    inx = ITIs2 > ITIMax  | ITIs2 < ITIMin;
    ITIs2(inx) = random('Normal',mng2,sdg,1,sum(inx));
end
ITIs3 = random('Uniform',ITIMin,ITIMax,1,NumTrials);
prr = rand(1,NumTrials);
rr = zeros(3,NumTrials);
rr(1,prr<pmx1) = 1;
rr(2,prr>=pmx1&prr<(pmx1+pmx2)) = 1;
rr(3,prr>=(pmx1+pmx2)) = 1;
ITIs = rr(1,:) .* ITIs1 + rr(2,:) .* ITIs2 + rr(3,:) .* ITIs3;
ITIMin = min(ITIs);
ITIMax = max(ITIs);

[nm xout] = hist(ITIs,cnts);
ivs = 10 ./ nm;
edges = [cnts-ivs; cnts+ivs];