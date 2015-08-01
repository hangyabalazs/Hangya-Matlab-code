function [H PerfGo PerfNoGo TimeGo TimeNoGo] = auditory_gonogo_psychplot3(mousename,session,lastn,valid_trials)
% AUDITORY_GONOGO_PSYCHPLOT2 plots psychometric performance plot for
% auditory go nogo task for multiple sessions.
% specify a folder with TE files.
% For all files in that folder,
% plot a psychometric, chronometric and d' curve.

% specify folder
if nargin < 1,
    mousename = 'mk05';
end

basedir = getpref('cellbase','datapath');

datapath =[basedir filesep mousename filesep session];

list=dir(datapath);
list=setdiff({list.name},{'.';'..';'.DS_Store'});
stc = strfind(list,'TE');
cf = ~cellfun(@isempty,stc);
list = list(cf);


if nargin < 3 
    list = list(1:end);
else
    if ~isempty(lastn)
        list = list(end-(lastn-1):end);
    else
        list = list(1:end);
    end
end
%
% list = listfiles(datapath);

% Specify colors.
colors = flipud(autumn(length(list)));
colors = flipud(jet(length(list)));
% colors = flipud(cool(length(list)));

% Specify plot properties.
ms = 8;
lw = 2;

H = figure;
set(gcf,'DefaultAxesFontName','Arial')
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'DefaultAxesBox','off')

% valid session ids.
valsess = [];

% For all files in that folder.
for iS = 1:length(list),
    TE = load([datapath filesep list{iS}]);
    if isequal(whichcb,'human_susattn')
        TE.StimulusDuration = log(TE.SoundIntensity);
    end
    
    if nargin == 4
        TE = filterTE(TE,valid_trials);
    end

    
    % Get session id.
    sessid{iS} = TE.sessionID{end-1};
    mouse_name{iS} = TE.Ratname{end-1};
    
    % Find unique Sound intensities used.
    [A, Ia, Ib] = unique(TE.StimulusDuration);
    
    NumStim = 0;
    if length(A) > NumStim,
        valsess = [valsess iS];
        % Initialization.
        PerfGo = nan(1,length(A));
        PerfNoGo = nan(1,length(A));
        TimeGo = nan(1,length(A));
        TimeNoGo = nan(1,length(A));
        
        % For each unique Sound Intensity.
        for i = 1:length(Ia),
            
            % Get performance on A+ trials.
            PerfGo(i) = nansum(TE.Hit(Ib==i))/(nansum(TE.Hit(Ib==i))+nansum(TE.Miss(Ib==i)));
            TimeGo(i) = nanmedian(TE.GoRT(Ib==i));
            
            % Performance on B- trials.
            PerfNoGo(i) = nansum(TE.FalseAlarm(Ib==i))/(nansum(TE.CorrectRejection(Ib==i))+nansum(TE.FalseAlarm(Ib==i)));
            %     NoGoPerf(i) = nansum(TE.CorrectRejection(Ib==i))/(nansum(TE.CorrectRejection(Ib==i))+nansum(TE.FalseAlarm(Ib==i)));
            TimeNoGo(i) = nanmedian(TE.NoGoRT(Ib==i));
            
        end
        
        subplot(251)
        plot(A,PerfGo,'LineStyle','-','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(iS,:),'Color',colors(iS,:),'MarkerSize',ms,'LineWidth',lw);
        hold on
        plot(A,PerfNoGo,'LineStyle','--','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(iS,:),'Color',colors(iS,:),'MarkerSize',ms,'LineWidth',lw);
        ylim([0 1])
        box off
        title('Performance')
        
        
        subplot(252)
        plot(A,TimeGo,'LineStyle','-','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(iS,:),'Color',colors(iS,:),'MarkerSize',ms,'LineWidth',lw);
        hold on
        plot(A,TimeNoGo,'LineStyle','--','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(iS,:),'Color',colors(iS,:),'MarkerSize',ms,'LineWidth',lw);
        box off
        title('Reaction time')
        xlabel('Signal intensity (dB)')
        
        SI = unique(TE.StimulusDuration);
        numbins = 100;
        EDGES = linspace(nanmin(TE.ReactionTime),nanmax(TE.ReactionTime),numbins);
        co=0;
        lw = 2;
        shadeness = linspace(0.9,0.1,length(SI));
        subplot(253)
        hold on
        for iSub = 1:length(SI),
            [N,X]= histc(TE.GoRT(TE.StimulusDuration==SI(iSub)),EDGES);
            h(iSub) = stairs(EDGES,cumsum(N)./max(cumsum(N)),'Color',[0 shadeness(iSub) 0],'LineWidth',lw);
            axis tight
        end
        
        % cum histogram of NoGoRT as function of SI
        subplot(254)
        hold on
        for iSub = 1:length(SI),
            [N,X]= histc(TE.NoGoRT(TE.StimulusDuration==SI(iSub)),EDGES);
            stairs(EDGES,cumsum(N)./max(cumsum(N)),'Color',[shadeness(iSub) 0 0],'LineWidth',lw);
            axis tight
        end
        
        subplot(255)
        hold on
        plot(A,PerfGo-PerfNoGo,'LineStyle','-','Marker','o','MarkerEdgeColor','none','MarkerFaceColor',colors(iS,:),'Color',colors(iS,:),'MarkerSize',ms,'LineWidth',lw);
        ylim([-1 1])
        title('d''')
        
        if iS == length(list),
%             legend(sessid(valsess),'Location','Best');
        end
    end
end

if length(unique(mouse_name))>1,
    warning('Data from more than one mouse')
    disp(unique(mouse_name))
else
    fstamp(unique(mouse_name))
end

TE = TrialEventsMerge(datapath,list);

bl = 18;
% matrix of blocks. No overlap.
mHit = reshape(padarray(TE.Hit,[0 bl-rem(length(TE.Hit),bl)],NaN,'post'),bl,[]);
mMiss = reshape(padarray(TE.Miss,[0 bl-rem(length(TE.Miss),bl)],NaN,'post'),bl,[]);
mCR = reshape(padarray(TE.CorrectRejection,[0 bl-rem(length(TE.CorrectRejection),bl)],NaN,'post'),bl,[]);
mFA = reshape(padarray(TE.FalseAlarm,[0 bl-rem(length(TE.FalseAlarm),bl)],NaN,'post'),bl,[]);

subplot(2,3,[4:6])
cla
% mGoPerf = nansum(mHit)./(nansum(mHit)+nansum(mMiss));
% mNoGoPerf = nansum(mFA)./(nansum(mFA)+nansum(mCR));
% plot(mGoPerf,'g.-');
% hold on
% plot(mNoGoPerf,'r.-');
% 
% plot((mGoPerf - mNoGoPerf),'k.-')
% 
% legend({'%Hit','%FA','d'''},'Location','best')
% % plot((mGoPerf - mNoGoPerf)./(mGoPerf + mNoGoPerf),'k.-')

[Bl, Ia, Ib] = unique(TE.BlockNum);
bGoPerf = nan(1,length(Bl));
for ibl = 1:length(Bl),
    bGoPerf(ibl) = nansum(TE.Hit(Ib==ibl))/(nansum(TE.Hit(Ib==ibl))+nansum(TE.Miss(Ib==ibl)));
    bNoGoPerf(ibl) = nansum(TE.FalseAlarm(Ib==ibl))/(nansum(TE.FalseAlarm(Ib==ibl))+nansum(TE.CorrectRejection(Ib==ibl)));
end

plot(bGoPerf,'go-');
hold on
plot(bNoGoPerf,'ro-');

plot((bGoPerf - bNoGoPerf),'ks-')

% lines for session numbers.
[Sess,Ie,junk] = unique(TE.SessionID);
line(repmat(TE.BlockNum(Ie),2,1),repmat(ylim',1,length(Ie)),'Color','k')
line(xlim,[0 0],'Color',[0.7 0.7 0.7])
yl = ylim;

% session ids. next to session starts.

[Sess,Is,junk] = unique(TE.SessionID,'first');
SessID = TE.sessionID(Ie); % because the first trial does not have a sessionID.

text(TE.BlockNum(Is),repmat(yl(1),length(Is),1),SessID,'Color','b','VerticalAlignment','top')

legend({'%Hit','%FA','d'''},'Location','best')

% RT vs ITI
figure
plot(TE.ITIDistribution,TE.GoRT,'.')
figure
SIn = SI(~isnan(SI));
maxint = SIn(end);
plot(TE.ITIDistribution(TE.StimulusDuration==maxint),TE.GoRT(TE.StimulusDuration==maxint),'.')
if length(SIn) > 1
    figure
    maxint2 = SIn(end-1);
    plot(TE.ITIDistribution(TE.StimulusDuration==maxint2),TE.GoRT(TE.StimulusDuration==maxint2),'.')
end

% Supress output
if nargout < 1
    clear
end

disp('hmm')



% -------------------------------------------------------------------------
function TE = filterTE(TE,valid_trials)

fnm = fieldnames(TE);
for k = 1:length(fnm)
    TE.(fnm{k}) = TE.(fnm{k})(valid_trials);
end