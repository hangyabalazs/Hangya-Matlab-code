function bimodalITI3
% 

% Specify folder
basedir = getpref('cellbase','datapath');
mousename = 'n029';
switch mousename
    case 'n023'
        
%         sessions = {'111220a' '111221a' '111222a' '111223a' '111227a' '111228a' ...
%             };    % bimodal, different mixing proportions (narrow peaks)
                
%         sessions = {'111229a' '111230a' '111231a' '120101a' '120102a' '120103a' ...
%             '120104a' '120105a' '120106a' '120110a' '120111a' '120112a' '120113a' ...
%             '120114a' '120115a' '120116a'};    % bimodal
        
%         sessions = {'120117a' '120118a' '120119a' '120120a' '120121a' '120122a'};   % unimodal
%         
%         sessions = {'120123a' '120124a' '120125a' '120126a' '120127a' '120128a'...
%             '120131a' '120201a' '120202a' '120203a'};   % unimodal with uniform
%         
        sessions = {'120207a' '120208a' '120209a' '120210a' '120211a' ...
            '120213a' '120214a' '120215a' '120216a' '120217a'};     % exponential
%         
%         sessions = {'120220a' '120221a'};       % unimodal again

    case 'n026'
        
        sessions = {'120207a' '120208a' '120209a' '120210a' '120211b' ...
            '120213a' '120214a' '120215a' '120216a'};     % exponential
        
%         sessions = {'120218a' '120220a' '120221a' '120222a' '120223a'};  % bimodal (not performing after Cosyne, not included) 
        
    case 'n027'
        
        sessions = {'120208a' '120209a' '120210a' '120211a'};     % exponential
        
%         sessions = {'120213b' '120214a' '120215a' '120215c' '120216a' '120217a' ...
%             '120218a' '120220a' '120221a' '120222a' '120223a' ...
%             '120301a' '120302a' '120303a' '120305a' '120306a'...
%             '120307a' '120308a' '120309a' '120310a' '120311a'...
%             '120313a' '120314a' '120315a'};     % bimodal (120215a and c: not doing)
        
    case 'n028'
        
        sessions = {'120207a' '120208a' '120210a' '120211a' ...
            '120213a' '120214a' '120215a' '120216a'};     % exponential
        
%         sessions = {'120217a' '120218a' '120220a' '120221a'...
%             '120222a' '120222c' '120223a' ...
%             '120301a' '120302a' '120303a' '120305a' '120306a' '120307a'};     % bimodal
        
    case 'n029'
        
        sessions = {'120204a' '120204b' '120206a' '120207a' '120207b' ...
            '120208a' '120209a' '120210a' '120211a' '120213a'...
            '120213b' '120214a' '120214b' '120215a' '120216b' '120217a' ...
            '120218a' '120220a' '120220b' '120221a' '120221b' '120222a' '120222b' ...
            '120301a' '120302a' '120303a' '120305a'...
            '120306a' '120307a' '120308a' '120309a' '120310a' '120311a'...
            '120313a' '120314a' '120315a'};     % bimodal (remark: 120209a session very long and good RT data!!!)
        
end

numSessions = length(sessions);

% Specify colors.
colors = flipud(autumn(length(numSessions)));
colors = flipud(jet(length(numSessions)));

% Specify plot properties.
ms = 8;
lw = 2;

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
    edges = 0.1:0.1:5;
%     edges = [0.1 0.3 0.8 1.5 1.9 2.1 3];
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    bno = length(cnts);
    maxint = SIn(end);
    if length(SI) > 1
        maxint2 = SIn(end-1);
    end
    [meanrt{iS} medianrt{iS} nrt{iS} allrt{iS} alliti{iS}] = dhist([TE.ITIDistribution' TE.GoRT'],edges);
    [meanhit{iS} medianhit{iS} nhit{iS} allhit{iS} alliti{iS}] = dhist([TE.ITIDistribution' TE.Hit'],edges);
    nhits{iS} = cellfun(@(g)sum(nan2zero(g)),allhit{iS});
    [meanmiss{iS} medianmiss{iS} nmiss{iS} allmiss{iS} alliti{iS}] = dhist([TE.ITIDistribution' TE.Miss'],edges);
    nmisss{iS} = cellfun(@(g)sum(nan2zero(g)),allmiss{iS});
    [~, ~, ~, allsi{iS}] = dhist([TE.ITIDistribution' TE.StimulusDuration'],edges);
    allrt2 = [allrt2 TE.GoRT];
    alliti2 = [alliti2 TE.ITIDistribution];
    allhit2 = [allhit2 TE.Hit];
    allmiss2 = [allmiss2 TE.Miss];
    allsi2 = [allsi2 TE.StimulusDuration];
%     [meanrt{iS} medianrt{iS} nrt{iS} allrt{iS} alliti{iS}] = ...
%         dhist([TE.ITIDistribution(TE.StimulusDuration==maxint)' TE.GoRT(TE.StimulusDuration==maxint)'],edges);
%     [meanrt{iS} medianrt{iS} nrt{iS} allrt{iS} alliti{iS}] = ...
%         dhist([TE.ITIDistribution(TE.StimulusDuration==maxint|TE.StimulusDuration==maxint2)'...
%         TE.GoRT(TE.StimulusDuration==maxint|TE.StimulusDuration==maxint2)'],edges);
    
%     figure
%     plot(TE.ITIDistribution,TE.GoRT,'.')
%     figure
%     
%     maxint = SIn(end);
%     plot(TE.ITIDistribution(TE.StimulusDuration==maxint),TE.GoRT(TE.StimulusDuration==maxint),'.')
%     if length(SI) > 1
%         figure
%         maxint2 = SIn(end-1);
%         plot(TE.ITIDistribution(TE.StimulusDuration==maxint2),TE.GoRT(TE.StimulusDuration==maxint2),'.')
%     end
end

mnsall = cell(1,length(cnts));
mnsmedian = nan(1,length(cnts));
for t = 1:length(cnts)
    for k = 1:numSessions
        mnsall{t} = [mnsall{t}; allrt{k}{t}];
        mnsall{t}(mnsall{t}<0.1) = NaN;
        mnsall{t}(mnsall{t}>0.45) = NaN;
    end
    mnsmedian(t) = nanmedian(mnsall{t});
end
figure
plot(cnts,mnsmedian)
hold on
[nm xout] = hist(TE.ITIDistribution,50);
stairs(xout,nm/300+0.2,'Color','r')

keyboard

nhi = cell2mat(nhits);
nmi = cell2mat(nmisss);
npe = nhi ./ (nhi + nmi);
nperf = sum(nhi) ./ (sum(nhi) + sum(nmi));
figure
plot(cnts,nperf)
hold on
[nm xout] = hist(TE.ITIDistribution,50);
stairs(xout,nm/200+0.5,'Color','r')

subst = 1:bno;
% subst = 12:29;
lensu = length(subst);
allrt_sorted = cell(lensu,1);
alln_sorted = [];
session_var = [];
cntr = 0;
for k = subst
    cntr = cntr + 1;
    for t = 1:numSessions
        atk = allrt{t}{k};
        atk = atk(~isnan(atk));
        allrt_sorted{cntr} = [allrt_sorted{cntr}; atk];
        session_var = [session_var; ones(length(atk),1)*t];
    end
    alln_sorted = [alln_sorted; ones(length(allrt_sorted{cntr}),1)*k];
end
figure
boxplot(cell2mat(allrt_sorted),alln_sorted)
anova1(cell2mat(allrt_sorted),alln_sorted)
rm_anova2(cell2mat(allrt_sorted),session_var,alln_sorted,alln_sorted(randperm(length(alln_sorted))),{'RT' 'RTT'})
lnt = length(alln_sorted);
rndrt = [ones(floor(lnt/2),1); 2*ones(ceil(lnt/2),1)];
rndrt = rndrt(randperm(lnt));
X = [cell2mat(allrt_sorted) alln_sorted rndrt session_var];
mixed_between_within_anova(X)
anova1(cell2mat(mnsall'),alln_sorted)

% linx = find((cnts>=0.125&cnts<0.525)|(cnts>1.525&cnts<2.675));
% sinx = find(~(cnts>=0.125&cnts<0.525)&~(cnts>1.525&cnts<2.675));
% linx = find((cnts<=0.65)|(cnts>=1.25&cnts<=1.85));
% sinx = find(~(cnts<=0.65)&~(cnts>=1.25&cnts<1.85));
linx = find((cnts<=0.55)|(cnts>=1.35&cnts<=1.9));
sinx = find((cnts>=0.65&cnts<=1.15)|(cnts>=2.45&cnts<2.9));
alliti_lin = [];
allhit_lin = [];
allmiss_lin = [];
allsi_lin = [];
allrt_lin = [];
for k = 1:numSessions
    alliti_lin = [alliti_lin; cell2mat(alliti{k}')];
    allhit_lin = [allhit_lin; cell2mat(allhit{k}')];
    allmiss_lin = [allmiss_lin; cell2mat(allmiss{k}')];
    allsi_lin = [allsi_lin; cell2mat(allsi{k}')];
    allrt_lin = [allrt_lin; cell2mat(allrt{k}')];
end
% linx2 = find((alliti_lin>=0.125&alliti_lin<0.525)|(alliti_lin>1.525&alliti_lin<2.675));
% sinx2 = find(~(alliti_lin>=0.125&alliti_lin<0.525)&~(alliti_lin>1.525&alliti_lin<2.675));
% linx2 = find((alliti_lin<=0.65)|(alliti_lin>=1.25&alliti_lin<=1.85));
% sinx2 = find(~(alliti_lin<=0.65)&~(alliti_lin>=1.25&alliti_lin<=1.85));
% linx2 = find((alliti_lin<=0.65)|(alliti_lin>=1.25&alliti_lin<=1.85));
% sinx2 = find((alliti_lin>=0.65&alliti_lin<=0.95)|(alliti_lin>=2.45&alliti_lin<=2.75));
% 
% linx2 = find((alliti2<=0.65)|(alliti2>=1.25&alliti2<=1.85));
% sinx2 = find((alliti2>=0.65&alliti2<=0.95)|(alliti2>=2.45&alliti2<=2.75));
linx2 = find((alliti2<=0.55)|(alliti2>=1.35&alliti2<=1.9));
sinx2 = find((alliti2>=0.65&alliti2<=1.15)|(alliti2>=2.45&alliti2<=2.9));

% linx2 = find(allrt2>nanmedian(allrt2));
% sinx2 = find(allrt2<nanmedian(allrt2));

% lhit = allhit_lin(linx2);
% lmiss = allmiss_lin(linx2);
% lsi = allsi_lin(linx2);
% lrt = allrt_lin(linx2);
% shit = allhit_lin(sinx2);
% smiss = allmiss_lin(sinx2);
% ssi = allsi_lin(sinx2);
% srt = allrt_lin(sinx2);
lhit = allhit2(linx2);
lmiss = allmiss2(linx2);
lsi = allsi2(linx2);
lrt = allrt2(linx2);
shit = allhit2(sinx2);
smiss = allmiss2(sinx2);
ssi = allsi2(sinx2);
srt = allrt2(sinx2);
numInt = length(SIn);
lPerfGo = nan(1,numInt);
lTimeGo = nan(1,numInt);
sPerfGo = nan(1,numInt);
sTimeGo = nan(1,numInt);
for psn = 1:numInt
    sn = SIn(psn);
    lPerfGo(psn) = nansum(lhit(lsi==sn)) / ...
        (nansum(lhit(lsi==sn)) + nansum(lmiss(lsi==sn)));
    lTimeGo(psn) = nanmedian(lrt(lsi==sn));
    
    sPerfGo(psn) = nansum(shit(ssi==sn)) / ...
        (nansum(shit(ssi==sn)) + nansum(smiss(ssi==sn)));
    sTimeGo(psn) = nanmedian(srt(ssi==sn));
end
figure
plot(lPerfGo,'Color',[0 0.8 0])
hold on
plot(sPerfGo,'Color',[0 1 0])
figure
plot(lTimeGo,'Color',[0 0.8 0])
hold on
plot(sTimeGo,'Color',[0 1 0])

ranksum(mnsmedian(linx),mnsmedian(sinx))    % sign
sn = 30;
l30 = lrt(lsi==sn);
s30 = srt(ssi==sn);
ranksum(l30(~isnan(l30)),s30(~isnan(s30)))  % sign for 11 sessions
sn = 40;
l40 = lrt(lsi==sn);
s40 = srt(ssi==sn);
ranksum(l40(~isnan(l40)),s40(~isnan(s40)))

rts1 = allrt2(linx2);
rts2 = allrt2(sinx2);
nanmedian(rts1)
nanmedian(rts2)
ranksum(rts1(~isnan(rts1)),rts2(~isnan(rts2)))
rperf1 = nansum(allhit2(linx2))/(nansum(allhit2(linx2))+nansum(allmiss2(linx2)))
rperf2 = nansum(allhit2(sinx2))/(nansum(allhit2(sinx2))+nansum(allmiss2(sinx2)))

keyboard

smns = smooth(mnsmedian,'linear',5);
figure;
plot(cnts,smns)
hold on
[nm xout] = hist(TE.ITIDistribution,50);
% stairs(xout,nm/300+0.2,'Color','r')
% stairs(xout,nm/900+0.18,'Color','r')
stairs(xout,nm/1000+0.25,'Color','r')

figure
hold on
for t = 1:length(alliti{k})
    for k = 1:numSessions
        plot(alliti{k}{t},allrt{k}{t},'.')
    end
end

figure
hold on
for k = 1:numSessions
    for t = 1:length(alliti{k})
        plot(alliti{k}{t},allrt{k}{t},'.')
    end
end

figure
hold on
for k = 1:numSessions
    plot(cnts,medianrt{k})
end

figure
mnn = cell2mat(medianrt);
plot(cnts,nanmean(mnn))

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