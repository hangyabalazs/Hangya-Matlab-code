%%

cellid = {'n071' '141218a'}
cellid = {'n038' '120922a'}
cellid = {'n013' '110526a'}
cellid = {'n023' '120106a'}
cellid = {'n046' '130108a'}
cellid = {'n072' '141222a'}
align_event = 'StimulusOn';
align_event = 'DeliverFeedback';
% align_event = 'LeftPortIn';

% Time window
wn = [-2 2];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% Calcualte PSTH
[psth, spsth, spsth_se, ~, spt] = ...
    ultimate_psth(cellid,'lick',align_event,wn,...
    'dt',dt,'sigma',0.02,'parts','#ResponseType','isadaptive',2,...
    'maxtrialno',Inf,'first_event','TrialStart','last_event','TrialEnd');

%%
% Lick raster
H = figure;
viewlick(cellid,'TriggerName',align_event,'SortEvent','TrialStart','eventtype','behav',...
    'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
    'Partitions','#ResponseType','window',[-5 5])
maximize_figure(H)

%%
% Replace PSTH
L = findobj(allchild(gca),'Type','line');   % lines in the legend
clr = flipud(get(L,'Color'));
cla
for k = 1:size(psth,1)
    errorshade(time,spsth(k,:),spsth_se(k,:),'LineWidth',2,...
        'LineColor',clr{k},'ShadeColor',clr{k})
    hold on
end
axis tight
set(gcf,'Renderer','OpenGL')

figure
errorshade(time,spsth(1,:),spsth_se(1,:),'LineWidth',2,...
    'LineColor',clr{1},'ShadeColor',clr{1})
errorshade(time,spsth(2,:),spsth_se(2,:),'LineWidth',2,...
    'LineColor',clr{2},'ShadeColor',clr{2})
hold on