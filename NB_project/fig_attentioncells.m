function fig_attentioncells
%FIG_ATTENTIONCELLS   Attention analysis figures.
%   FIG_ATTENTIONCELLS makes figures for reaction time prediction of
%   example neurons (PSTH, ROC; ROC window: 50 ms, causal).
%
%   See also REGRESSION_ANALYSIS, STIMRASTER2, NBITIFR and
%   COND_ACCURACY_FR3.

% PSTH, raster plot, ROC analysis
% main({'n023_111214b_7.1'},[],false)
% main({'n040_121110a_1.2'},[],false)   % extended data fig example 1
% main({'n028_120217a_7.1'},[],false)  % main fig.2
% main({'n028_120302a_7.1'},[],false)   % extended data fig example 2
main({'n020_111021a_4.1'},[],false)   % extended data fig example 3

% -------------------------------------------------------------------------
function [error_list1 error_list2] = main(cellids,resdir,issave)
% Based on NBSTIMRASTER2 and NBITIFR

% PSTH aligned to trial start
% [allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
%     error_list1] = poppsth(cellids,'ITI',[],'StimulusOn',resdir,issave); %#ok<*ASGLU>
% time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
% fnm = [resdir 'allPSTH_ITI.mat'];
% if issave
%     save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
%         'tags','error_list1')
% end

% PSTH aligned to stimulus on
[allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
    error_list2] = poppsth(cellids,'StimulusOn','LastITIBegins',[],resdir,issave); %#ok<*ASGLU>
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
fnm = [resdir 'allPSTH_StimulusOn.mat'];
if issave
    save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
        'tags','error_list2')
end

% -------------------------------------------------------------------------
function [allpsth allpsth2 allspsth allspsth2 allstats ...
    ok_ids wn error_list] = poppsth(cellids,alignevent,...
    first_event,last_event,resdir,issave)

% Time window
wn = [-2.5 1];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector
clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];

% Call 'main'
% Preallocate
[allpsth allspsth] = deal([]);
allstats = struct([]);
error_list = struct('cellid',{},'message',{},'exception',{});  % keep a list of caught errors 
errinx = 0;
ok_ids = {};
NumCells = length(cellids);
for k = 1:NumCells
    cellid = cellids{k};   % cell ID
    disp(cellid)
    
    try
        
        % Check for 'DeliverFeedback' event
%         alignfilter = 'Hit==1';
%         alignevent = 'ITIback';   % align to stimulus on but goes back only to trial start
        
        % Calcualte PSTH
        [psth_shortRT, spsth_shortRT, spsth_se_shortRT, ~, spt_shortRT, stats_shortRT] = ...
            ultimate_psth(cellid,'trial',alignevent,wn,...
            'dt',dt,'display',false,'sigma',0.05,'parts','all','isadaptive',2,...
            'event_filter','selectGoRT','filterinput',[0.1 0.5],'maxtrialno',Inf,...
            'first_event',first_event,'last_event',last_event,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        [psth_longRT, spsth_longRT, spsth_se_longRT, ~, spt_longRT, stats_longRT] = ...
            ultimate_psth(cellid,'trial',alignevent,wn,...
            'dt',dt,'display',false,'sigma',0.05,'parts','all','isadaptive',2,...
            'event_filter','selectGoRT','filterinput',[0.5 0.9],'maxtrialno',Inf,...
            'first_event',first_event,'last_event',last_event,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        figure
        errorshade(time,spsth_shortRT,spsth_se_shortRT,'LineWidth',2,...
            'LineColor',clr(1,:),'ShadeColor',clr(1,:))
        hold on
        errorshade(time,spsth_longRT,spsth_se_longRT,'LineWidth',2,...
            'LineColor',clr(4,:),'ShadeColor',clr(4,:));
        set(gca,'Color',[0.8 0.8 0.8])
        
        % Concatenate data from different cells
        psth = [psth_shortRT; psth_longRT];
        spsth = [spsth_shortRT; spsth_longRT];
        stats = [stats_shortRT; stats_longRT];
        allpsth = [allpsth; psth]; %#ok<AGROW>
        allspsth = [allspsth; spsth]; %#ok<AGROW>
        allstats = [allstats; stats]; %#ok<AGROW>
        ok_ids{end+1} = cellid; %#ok<AGROW>
        
        % ROC analysis
        wns = 0.05 / dt;   % 50 ms window (in data points)
        sht = 0.02 / dt;   % shift between overlapping windows (in data points)
        nmw = floor((diff(wn)/dt-wns)/sht) + 1;   % number of overlapping windows
        [ROC P SE] = deal(nan(1,nmw));
        for w = 1:nmw
            inx1 = (w - 1) * sht + 1;
            inx2 = inx1 + wns;
            [ROC(w) P(w) SE(w)] = rocarea(nansum(spt_shortRT(:,inx1:inx2),2),...
                nansum(spt_longRT(:,inx1:inx2),2),'transform','scale','bootstrap',1000);
        end
        ROCtime = (wn(1)/dt+wns/2:sht:wn(2)/dt-wns/2) * dt;   % time vector
        ROCtime = ROCtime + wns * dt / 2;   % make the window causal
        ROC = smooth(nan2zero(ROC),'linear',5);   % NaN when all spike counts are 0 for both distributions
        
        figure   % plot
%         plot(ROCtime,smooth(ROC,'linear',11),'k')
        errorshade(ROCtime(ROCtime>=-2&ROCtime<=0.6),ROC(ROCtime>=-2&ROCtime<=0.6),SE(ROCtime>=-2&ROCtime<=0.6),...
            'LineColor','k','ShadeColor','k')
        
        % Save figure
        Hroc = gcf;
        if issave
            fnm = [resdir cellidt '_' alignevent '_ROC.fig'];
            saveas(Hroc,fnm)
            fnm = [resdir cellidt '_' alignevent '_ROC.jpg'];
            saveas(Hroc,fnm)
            fnm = [resdir cellidt '_' alignevent '_ROC.mat'];
            save(fnm,'ROCtime','ROC')
        end
%         close(Hroc)
        
        % Raster plot with sorted trials
        Hr2 = figure;
        viewcell2b(cellids(k),'TriggerName',alignevent,'SortEvent','LeftPortIn',...
            'eventtype','behav','ShowEvents',{{'LeftPortIn'}},...
            'Partitions','#StimulusID&Hit','Num2Plot',100,'window',wn)
        maximize_figure(Hr2)
%         ln = findobj(allchild(gca),'type','line','color','k');
%         set(ln,'LineWidth',1)
%         setmyplot_balazs
%         set(gco,'FaceColor',[0.9490 0.8000 0.0471],'EdgeColor',[0.9490 0.8000 0.0471])
%         set(gco,'Color',[1.0000    0.6941    0.3922],'LineWidth',2)
        cla   % replace PSTH
        darkgreen = [0 0.5 0];
        lightgreen = [0 1 0];
        errorshade(time,spsth_shortRT,spsth_se_shortRT,'LineWidth',2,...
            'LineColor',darkgreen,'ShadeColor',darkgreen)
        hold on
        errorshade(time,spsth_longRT,spsth_se_longRT,'LineWidth',2,...
            'LineColor',lightgreen,'ShadeColor',lightgreen);
        set(gca,'Color',[0.8 0.8 0.8])
        axis tight
        line([0 0],ylim,'color',[0.7412 0.0078 0.0078],'LineWidth',1)
        legend({'short RT','long RT'})
        set(gcf,'Renderer','OpenGL')
        A = findobj(allchild(gcf),'Type','axes');  % overlay zero line
        axes(A(5));
        line([0 0],ylim,'LineWidth',3,'Color',[0 0 0.75])
        fnm = [resdir cellidt '_' alignevent '_rasterPSTH.jpg'];   % save
%         saveas(Hr2,fnm)
        if issave
            set(Hr2,'PaperPositionMode','auto')
            set(Hr2,'InvertHardcopy','off')
            print(Hr2,'-djpeg',fnm)
        end
%         close(Hr2)
        
    catch ME
        
        % Error handling
        errinx = errinx + 1;  % error counter
        error_list(errinx).cellid = cellid;   % record problematic cell ID
        error_list(errinx).message = ME.message;   % error message
        error_list(errinx).exception = ME;   % exception structure
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
    end
end

% Standardize
allpsth2 = allpsth;
allspsth2 = allspsth;
for k = 1:size(allpsth,1)
    allpsth2(k,:) = standardize(allpsth(k,:));
    allspsth2(k,:) = standardize(allspsth(k,:));
end