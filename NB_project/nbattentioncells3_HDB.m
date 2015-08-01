function nbattentioncells3
%NBATTENTIONCELLS3   Attention analysis.
%   NBATTENTIONCELLS3 performs various analyses to asess whether cells are
%   related to accuracy, action (animal's response), foreperiod
%   distribution or reaction time. See details about the analyses in
%   REGRESSION_ANALYSIS, NBSTIMRASTER2, NBITIFR and COND_ACCURACY_FR3.
%
%   See also REGRESSION_ANALYSIS, NBSTIMRASTER2, NBITIFR and
%   COND_ACCURACY_FR3.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   19-Dec-2013

%   Edit log: BH 12/19/13

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m
% load('C:\Balazs\_analysis\NB\regression\All\new2\regression_results.mat')

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m; trials from the lower 15 percentile of the RT
% distribution excluded
% load('C:\Balazs\_analysis\NB\regression\All\impexclude\regression_results.mat')

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m; trials from the lower and upper 10 percentile
% of the RT distribution excluded
global DATAPATH
% load([DATAPATH 'NB\regression\All\slowfastexclude\regression_results.mat'])
load([DATAPATH 'HDB\regression_newdata\All\slowfastexclude\regression_results.mat'])

% Directories
% resdir = fullfile(DATAPATH,'NB','attentioncells3_50ms_',filesep);   % results directory
resdir = fullfile(DATAPATH,'HDB','attentioncells3_newdata',filesep);   % results directory
issave = true;

% p-values
NumCells = length(p);
[p2 R2] = deal(nan(1,NumCells));
for k = 1:NumCells
    if ~isempty(p(k).iti05)
        p2(k) = p(k).iti05;
        R2(k) = R(k).iti05;
    end
end
p = p2;
R = R2;

% Areas
HDB = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);
[~, HDBinx] = intersect(I,HDB);

% Cell types
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
[~, ChATinx] = intersect(I,ChAT);

% NB cells with significant correlation (p<0.01)
nninx = find(~isnan(p));
HDBbinx = intersect(HDBinx,nninx);
HDBsinx = intersect(HDBinx,find(p<0.01));
CORR = I(HDBsinx);

% Percentage of significantly activated cells for different areas
pct.HDB = sum(p(HDBinx)<0.01) / sum(~isnan(p(HDBinx)));

% Percentage of significantly activated cells for different cell types
pct.ChAT = sum(p(ChATinx)<0.01) / sum(~isnan(p(ChATinx)));

% Save
if issave
    fnm = [resdir 'pct.mat'];
    save(fnm,'p','R','pct')
end

% PSTH, raster plot, ROC analysis
main(I(HDBsinx),resdir,issave)

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

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
% list_name1 = ['error_list1_' dsr];  % variable name
% list_name1 = regexprep(list_name1,'-','_');
% list_name1 = regexprep(list_name1,' ','_');
% assignin('base',list_name1,error_list1)   % assign error list in base workspace
% error_fname1 = fullfile(resdir,[list_name1 '.mat']);   % file name
% save(error_fname1,'error_list1')   % save error list
list_name2 = ['error_list2_' dsr];  % variable name
list_name2 = regexprep(list_name2,'-','_');
list_name2 = regexprep(list_name2,' ','_');
assignin('base',list_name2,error_list2)   % assign error list in base workspace
error_fname2 = fullfile(resdir,[list_name2 '.mat']);   % file name
save(error_fname2,'error_list2')   % save error list

% -------------------------------------------------------------------------
function [allpsth allpsth2 allspsth allspsth2 allstats ...
    ok_ids wn error_list] = poppsth(cellids,alignevent,...
    first_event,last_event,resdir,issave)

% Time window
wn = [-3 3];   % in seconds
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
            'dt',dt,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','selectGoRT','filterinput',[0.1 0.5],'maxtrialno',Inf,...
            'first_event',first_event,'last_event',last_event,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        [psth_longRT, spsth_longRT, spsth_se_longRT, ~, spt_longRT, stats_longRT] = ...
            ultimate_psth(cellid,'trial',alignevent,wn,...
            'dt',dt,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','selectGoRT','filterinput',[0.5 0.9],'maxtrialno',Inf,...
            'first_event',first_event,'last_event',last_event,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
%         figure
%         errorshade(time,spsth_shortRT,spsth_se_shortRT,'LineWidth',2,...
%             'LineColor',clr(1,:),'ShadeColor',clr(1,:))
%         hold on
%         errorshade(time,spsth_longRT,spsth_se_longRT,'LineWidth',2,...
%             'LineColor',clr(4,:),'ShadeColor',clr(4,:));
%         set(gca,'Color',[0.8 0.8 0.8])
        
        % Concatenate data from different cells
        psth = [psth_shortRT; psth_longRT];
        spsth = [spsth_shortRT; spsth_longRT];
        stats = [stats_shortRT; stats_longRT];
        allpsth = [allpsth; psth]; %#ok<AGROW>
        allspsth = [allspsth; spsth]; %#ok<AGROW>
        allstats = [allstats; stats]; %#ok<AGROW>
        ok_ids{end+1} = cellid; %#ok<AGROW>
        
        % Save figure
%         H = gcf;
        cellidt = regexprep(cellid,'\.','_');
%         if issave
%             fnm = [resdir cellidt '_' alignevent '_PSTH.fig'];
%             saveas(H,fnm)
%         end
%         close(H)
        
        % ROC analysis
%         wns = 0.2 / dt;   % 200 ms window (in data points)
        wns = 0.05 / dt;   % 50 ms window (in data points)
        sht = 0.02 / dt;   % shift between overlapping windows (in data points)
        nmw = floor((diff(wn)/dt-wns)/sht) + 1;   % number of overlapping windows
        ROC = nan(1,nmw);
        for w = 1:nmw
            inx1 = (w - 1) * sht + 1;
            inx2 = inx1 + wns;
            ROC(w) = rocarea(nansum(spt_shortRT(:,inx1:inx2),2),...
                nansum(spt_longRT(:,inx1:inx2),2),'transform','scale');
        end
        ROCtime = (wn(1)/dt+wns/2:sht:wn(2)/dt-wns/2) * dt;   % time vector (made causal in versions of fig_attentioncells)
        
        figure   % plot
        plot(ROCtime,ROC)
        
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
        close(Hroc)
        
        % Raster plot with sorted trials
        Hr2 = figure;
        viewcell2b(cellids(k),'TriggerName',alignevent,'SortEvent','StimulusOff',...
            'eventtype','behav','ShowEvents',{{'StimulusOn','StimulusOff','LastITIBegins'}},...
            'Partitions','#StimulusID&Hit','Num2Plot',100,'window',wn)
        maximize_figure(Hr2)
        cla   % replace PSTH
        errorshade(time,spsth_shortRT,spsth_se_shortRT,'LineWidth',2,...
            'LineColor',clr(1,:),'ShadeColor',clr(1,:))
        hold on
        errorshade(time,spsth_longRT,spsth_se_longRT,'LineWidth',2,...
            'LineColor',clr(4,:),'ShadeColor',clr(4,:));
        set(gca,'Color',[0.8 0.8 0.8])
        axis tight
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
        close(Hr2)
        
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