function auditory_gonogo_perf_rt(animalname,Numsess)
% auditory_gonogo
% SPR 2011/02/14
% create psychometric performance plot with mean of multiple sessions
% does not separate light stimulation trials
% create a cell array with each analyzed session file name
% startpath='C:\ratter\SoloData\Data\hanna';

basedir = '/Users/ranades/Documents/work/Data/behavior/auditory_gonogo';

if nargin<1,
    animalname = input('Enter the animal name for the session(s) you wish to analyze\n','s');
end

switch animalname(1:2)
    case 'jw'
%         datapath=['/Users/ranades/SoloData/Data/jordan' filesep animalname];
        datapath=[ basedir filesep animalname];
    case 'mk'
%         datapath=['/Users/ranades/SoloData/Data/matt'filesep animalname];
        datapath=[ basedir filesep animalname];
    case 'nb'
%         datapath=['/Users/ranades/SoloData/Data/matt'filesep animalname];
        datapath=[ basedir filesep animalname];
    otherwise
        error('unknown animal name')
end

% specify which sessions to plot.
if nargin<1,
    try
        file_list=uigetfile('TE*.mat','Select extracted sessions to analyze',[startpath filesep animalname filesep],'MultiSelect','on');
    catch
    %     errordlg('You have not entered a valid animal name. Please select an animal from the directory.');
    %     animal_directory = uigetdir('C:\ratter\SoloData\Data\hanna','Select animal');
    %     animal_name = char(animal_directory(31:length(animal_directory)));
    %     file_list=uigetfile(animal_directory,'Select extracted sessions to
    %     analyze','MultiSelect','on');
    end
    if~iscell(file_list) % ensures that file_list is a cell array
        file_list={file_list};
    end
else
    % list of all trial events files.
    list=dir(datapath);
    list=setdiff({list.name},{'.';'..';'.DS_Store'});
    
    if isnumeric(Numsess),
        file_list = list(end-(lastn-1):end);
    elseif ischar(Numsess),
        switch Numsess
            case 'last'
                
            case 'all'
                file_list = list;
                
            case 'YFF' % your favorite session filter.
                
            otherwise
                error(['Unkown session filter type' Numsess])
        end
    else
        error(['Unknown data type' Numsess])
    end % session specifier 
end % nargin


% create a loop that calculates PsychPerf at each stimdur for each session

NumSess = length(file_list);

StimDurList=[];
g.plotearlyresp=0;
g.plotlightstim=1;

for i=(1:NumSess)
    %determine filepath and load data
    filename=char(file_list(i));
    filepath = [datapath filesep filename];
    load(filepath)
    GoPerf(i) = nansum(Hit)/(nansum(Hit)+nansum(Miss));
    NoGoPerf(i) = nansum(FalseAlarm)/(nansum(FalseAlarm)+nansum(CorrectRejection));
    GoRT_distrib{i} = ReactionTime(Hit == 1);
    NoGoRT_distrib{i} = ReactionTime(FalseAlarm == 1);
    NumRewards(i) = nansum(Hit);
    NumResponses(i) = nansum(FalseAlarm) + nansum(Hit);
end

%%
figure(1)

subplot(221)
cla
plot(GoPerf,'go-')
hold on
plot(NoGoPerf,'ro-')
ylim([0 1])
title('learning')

subplot(222)
cla
plot(GoPerf-NoGoPerf,'bo-')
ylim([0 1])
title('learning d''')

subplot(223)
cla
GoRTstats = nans(NumSess,2);
NoGoRTstats = nans(NumSess,2);
for iS = 1:NumSess,
    GoRTstats(iS,:) = [nanmean(GoRT_distrib{iS}) nanstd(GoRT_distrib{iS}) ];
    NoGoRTstats(iS,:) = [nanmean(NoGoRT_distrib{iS}) nanstd(NoGoRT_distrib{iS})];
end
errorbar(1:size(GoRTstats,1),GoRTstats(:,1),GoRTstats(:,2),'go-')
hold on
errorbar(1:size(NoGoRTstats,1),NoGoRTstats(:,1),NoGoRTstats(:,2),'ro-')
title('RT')

subplot(224)
title('motivation')
cla
plot(1:NumSess,NumRewards,'bo-')
hold on
plot(1:NumSess,NumResponses,'ko-')
% legend('Rewards','tries')
text(0.95,0.95,num2str(NumRewards(end)),'Units','normalized','VerticalAlignment','top','FontSize',16,'HorizontalAlignment','right')
fstamp(animalname,[],'top-center')


%%
if 0,
    
Perf = cell(NumSess,1);
NoGoPerf = cell(NumSess,1);
StimList = cell(NumSess,1);
figure(3)
cla
hold on
for i=(1:NumSess)
    %determine filepath and load data
    filename=char(file_list(i));
    filepath = [datapath filesep filename];
    load(filepath)
    stims = unique(StimulusDuration);
    for iS = 1:length(stims),
        if stims(iS) ~= 0,
            Perf{i} = [Perf{i} nansum(Hit(StimulusDuration==stims(iS)))/...
                (nansum(Hit(StimulusDuration==stims(iS)))+...
                nansum(Miss(StimulusDuration==stims(iS))))];
        else
            Perf{i} = [Perf{i} nansum(FalseAlarm(StimulusDuration==stims(iS)))/...
                (nansum(FalseAlarm(StimulusDuration==stims(iS)))+...
                nansum(CorrectRejection(StimulusDuration==stims(iS))))];
        end
        StimList{i} = [StimList{i} stims(iS)]; 
    end
    disp(Perf{i})
%     plot(randn+StimList{i}(2:end),Perf{i}(2:end),'ro-')
    plot(StimList{i}(2:end),Perf{i}(2:end),'ko-')
%                 
%     GoRT_distrib{i} = ReactionTime(Hit == 1);
%     NoGoRT_distrib{i} = ReactionTime(FalseAlarm == 1);
end
    ylim([0 1])
end
% disp(Perf)
%%
% 
% 
% figure(1)
% 
% subplot(221)
% cla
% plot(GoPerf,'go-')
% hold on
% plot(NoGoPerf,'ro-')
% ylim([0 1])
% title('learning')
% 
% subplot(222)
% cla
% plot(GoPerf-NoGoPerf,'bo-')
% ylim([0 1])
% title('learning d''')
% 
% subplot(223)
% cla
% GoRTstats = nans(NumSess,2);
% NoGoRTstats = nans(NumSess,2);
% for iS = 1:NumSess,
%     GoRTstats(iS,:) = [nanmean(GoRT_distrib{iS}) nanstd(GoRT_distrib{iS}) ];
%     NoGoRTstats(iS,:) = [nanmean(NoGoRT_distrib{iS}) nanstd(NoGoRT_distrib{iS})];
% end
% errorbar(1:size(GoRTstats,1),GoRTstats(:,1),GoRTstats(:,2),'go-')
% hold on
% errorbar(1:size(NoGoRTstats,1),NoGoRTstats(:,1),NoGoRTstats(:,2),'ro-')
% title('RT')
% 
% fstamp(animalname,[],'top-center')
% 
