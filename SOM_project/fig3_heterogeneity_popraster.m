function [ R2 ] = fig3_heterogeneity_popraster(cellids)
%FIG3_HETEROGENEITY_POPRASTER get overall response profile for various
%events.
%   This analysis will use a event-spike cross correlation measure to
%   calculate selectivity of the neuron to various events in the task.
%   Event aligned binraster will be used for the cross correlation, and
%   compared with the shuffle corrected version to determine significance.

if nargin < 1,
    load('tagged_list_agreed');
    cellids = {'d029_100224b_1.1' , 'd042_100722a_6.1' , 'd042_100729a_1.3' , 'd043_100721a_3.3',  'd043_100721a_5.1', 'd043_100722a_5.1', 'd043_100731a_6.1' ,...
    'd043_100723a_2.2'    'd043_100731a_1.1'    'd043_100722a_1.1'    'd046_100816a_3.1'    'd046_100816a_5.2'    'd046_100818a_5.6',...
    'd035_100430a_1.2'    'd042_100803a_1.2'    'd047_100812a_1.1'    'd047_100812a_1.2',...
    'd042_100729a_3.1'    'd042_100729a_6.1'    'd043_100801a_4.1'    'd043_100802a_6.1',...
    'd042_100722a_3.1', 'd042_100809b_3.1' , 'd042_100809b_5.1', 'd043_100720a_2.1', 'd043_100722a_6.3' , 'd043_100723a_3.2','d047_100815a_2.1',...
    'd047_100808a_6.1'    'd047_100808a_6.4'    'd043_100731a_5.1'};
%     cellids = pv_beh;
    

% % %     
% %     
%     cellids = som_beh([1:19 21:end]);
end

% 'TriggerZoneIn','TriggerZoneOut','ReminderCue','Zone1FirstExit','HomeZoneIn1SpeedB1','HomeZOneOut1'
% decide the trigger events, and previous and nextevent. 
% ONLY FIXED EVENTS!
% iE = 1;
% events{iE} = {'TriggerZoneIn' 'PreviousHomeZoneOut' 'TriggerZoneOut'};
% iE = iE + 1;
% events{iE} = {'TriggerZoneOut' 'TriggerZoneIn' 'HomeZoneIn1'};
% iE = iE + 1;
% events{iE} = {'ReminderCue' 'HomeZoneIn1' 'HomeZoneOut1'};
% iE = iE + 1;
% events{iE} = {'HomeZoneIn1' 'TriggerZoneOut' 'HomeZoneOut1'};
% iE = iE + 1;
% events{iE} = {'WaterValveOn' 'TriggerZoneOut' 'HomeZoneOut1'};
% iE = iE + 1;
% events{iE} = {'HomeZoneOut1' 'HomeZoneIn1' 'NextTriggerZoneIn'};

iE = 1;
events{iE} = {'TriggerZoneIn'};
iE = iE + 1;
events{iE} = {'TriggerZoneOut'};
iE = iE + 1;
events{iE} = {'WaterValveOn'};
iE = iE + 1;
% events{iE} = {'ReminderCue'};
% iE = iE + 1;
% events{iE} = {'Zone1FirstExit'};
% iE = iE + 1;
events{iE} = {'HomeZoneIn1SpeedB1'};
iE = iE + 1;
events{iE} = {'HomeZoneOut1'};
% 

win = [-0.5 0.5];
g.eventtype = 'behav';
g.sigma = 0.04;
g.window = win;
g.Partitions = 'all';
win_margin = [0 0];
g.dt = 0.02;

margin = g.sigma*3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array

NumCells = length(cellids);
NumEvents = length(events);
% setmyfigure(gcf)
% sub_h = set_subplots(1,NumEvents,0,0);

% figure(39)
% clf
% setmyfigure(gcf)
% sub_h = set_subplots(2,ceil(length(events)/2),0.1,0.1);
% for each cell,
for iCell = 1:NumCells,
    cellid = cellids(iCell);
    
    % load prealigned spikes.
    switch g.eventtype
    case 'stim'
        TE=loadcb(cellid,'StimEvents');
        SP=loadcb(cellid,'STIMSPIKES');
    case {'event','behav'}
        TE=loadcb(cellid,'TrialEvents');
        SP=loadcb(cellid,'EVENTSPIKES');    
    end
    % for each event,
    for iE = 1:NumEvents,
        
        TriggerEvent = events{iE}{1};
%         PrevEvent = events{iE}{2};
%         NextEvent = events{iE}{3};
%         
        trigger_pos = findcellstr(SP.events(:,1),TriggerEvent);
        
%         % get valid windows for Previous and Next Events.
%         prev_windows = get_last_evtime(TE,TriggerEvent,PrevEvent);
%         next_windows = get_last_evtime(TE,TriggerEvent,NextEvent);
%         ev_windows = [prev_windows(:,2) next_windows(:,2)];
%         
%         % get valid trials based on windows.
%         % wierd isnan accounts for either window being a NaN
%         % if Previous event came after triggerevent and 
%         % if next event came before triggerevent.        
%         valid_i = sum(~isnan(ev_windows),2)==2 & ...
%                     ev_windows(:,1) < 0 & ...
%                     ev_windows(:,2) > 0;
        valid_i = 1:length(SP.event_stimes{trigger_pos});
        
        stimes  = SP.event_stimes{trigger_pos}(valid_i);
%         ev_windows = ev_windows(valid_i,:);
        
        % make the raster.
%         binraster = stimes2binraster(stimes,time,g.dt,ev_windows,win_margin);
        
        binraster = stimes2binraster(stimes,time,g.dt);
        
        [COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
        
        [psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_i);
        
        
        timearray_inx = time>g.window(1) & time < g.window(2);
        time1 = time(timearray_inx);
        spsth = spsth((timearray_inx));
                
        RESPMATRIX(:,iE,iCell) = spsth;

%         imagesc(time,1:sum(valid_i),binraster,'Parent',sub_2(1));
%         colormap(flipud(gray))
%         title(sub_2(1),TriggerEvent)
%         axes(sub_2(2));
%         stairs(time,nanmean(binraster,1))
%         set(sub_2(1),'Tag','raster')
%         set(sub_2(2),'Tag','psth')
        
    end
    % normalize by zscore
    
    % plot as a heat map.
    
%     ylims = cell2mat(get(findobj('Tag','psth'),'YLim'));
%     new_ylim = [min(ylims(:,1)) max(ylims(:,2))];
%     set(findobj('Tag','psth'),'YLim',new_ylim)
%     fstamp(cellid,10,'bottom-right')
%     pause(5)
end

R2 = reshape(RESPMATRIX,numel(RESPMATRIX)/NumCells,NumCells)';
for k = 1:size(R2,1)
    R2(k,:) = (R2(k,:) - min(R2(k,:))) / (max(R2(k,:)) - min(R2(k,:)));
end
% R2 = zscore(R2)';

figure
set(gcf,'Position',[360   274   815   404]);
clf
sub_h = set_subplots(1,NumEvents,0.01,0.1);
start_inx = 1:length(time1):size(R2,2);
for iP = 1:length(sub_h),
    subplot(sub_h(iP))
    imagesc(R2(:,start_inx(iP):length(time1)+start_inx(iP)-1));
    set(gca,'CLim',[0 1])
    set(gca,'XGrid','Off','YGrid','On','YTick',[0.5:NumCells+0.5],...
        'YTickLabel',[],'XTick',[1 ceil(length(time1)/2) length(time1)],'XTickLabel',sort([g.window 0]),'TickDir','out')
    set(gca,'TickLength',[0 0])
    line([find(time1==0,1,'first') find(time1==0,1,'first')],ylim,'Color','k')
    xlabel(events{iP})
%     axis off

end

% function binraster = get_binraster(Trig,Prev,Next,win,dt);
%  GET_BINRASTER for given trigger event and previous and next event, 
    




% % TriggerName mismatch
% if (trigger_pos == 0)
%     error('Trigger name not found');
% else
%     g.TriggerEvent=SP.events{trigger_pos,2};
% end
% 
% if ~isfield(TE,g.TriggerEvent),
%     error('TriggerEvent mismatch: supply correct Events structure')
% end
% 
% 
% g.TriggerEvent=SP.events(trigger_pos,2);

% if ~iscellstr(g.LastEvents) & (strcmpi(g.LastEvents,'none') | isempty(g.LastEvents))
%     window_margin = SP.events{trigger_pos,4};
%     ev_windows = SP.event_windows{trigger_pos};
% else
%     window_margin = [g.window(1)-2*g.dt 0];
%     ev_windows = get_last_evtime(TE,g.TriggerEvent,g.LastEvents);
% end
% 
% 
% alltrials = 1:size(SP.event_stimes{1},2);
% stimes  = SP.event_stimes{trigger_pos}(alltrials);
% windows = SP.event_windows{trigger_pos}(:,alltrials);
% 
% % restrict spikes from flanking events.
%     switch g.spikerestrict
%         case 'none'
%             
%         case 'both'
%             windows = [TE.(g.PrevEvent) - TE.(g.TriggerEvent); TE.(g.NextEvent) - TE.(g.TriggerEvent)];
%         
%         case 'pre'
%             windows(1,:) = TE.(g.PrevEvent) - TE.(g.TriggerEvent);
%         
%         case 'post'
%             windows(2,:) = TE.(g.NextEvent) - TE.(g.TriggerEvent);
%     end
% %     margins = [0 0];
%     
%  
% ev_windows = windows;
% 
% % margin = 	*3;     % add an extra margin to the windows
% % time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array
% 
% %%% MAKE THE MAIN RASTER
% binraster = stimes2binraster(stimes,time,g.dt,ev_windows,window_margin);
% 

end