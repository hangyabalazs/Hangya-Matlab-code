function error_list = nbstimraster(I,issave)
%NBSTIMRASTER   Peri-event time histogram.
%   NBSTIMRASTER(I,ISSAVE) calculates 'doubly adaptive' PSTHs for a
%   set of cells (see DAPSTH) aligned to 'Go' or 'No-go' response onset
%   ('Hit' or 'False Alarm'). Raster plots are also generated. 
%   Input parameters:
%       I - list of cell IDs or index set to CELLIDLIST  (see CellBase 
%           documentation); if empty or not specified, all well-separated
%           cells are selected (ID>20, L-ratio<0.15; see LRATIO) from basal
%           forebrain areas
%       ISSAVE - controls saving
%
%   ERROR_LIST = NBSTIMRASTER(I,ISSAVE) returns a structure with all
%   caught errors. ERROR_LIST includes cell IDs, error messages and the
%   captured exception variables. A time stamped ERROR_LIST is also
%   assigned in base workspace and saved to the results directory
%   automatically.
%
%   See also NBRESPNSESORTER, DAPSTH, ULTIMATE_PSTH, PSTH_STATS, LRATIO,
%   GAP_STATISTICS and NBPOPPSTH_CALL.

%   Edit log: BH 7/3/12, 8/13/12, 4/23/12
 
% Pass the control to the user in case of error
dbstop if error
if ~isequal(whichcb,'NB')
    choosecb('NB')
end

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = false;
end
if nargin < 1
    I = [];
end
 
% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'attentioncells' fs];   % results directory

% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');
 
% List of cellIDs
if isempty(I)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    vldty = getvalue('validity');
    area1 = getvalue('Area1');
    area2 = getvalue('Area2');
    inb = isnb(area1,area2);   % select NB cells
    ptinx = vldty == 1 & ID > 20 & Lratio < 0.15 & inb;   % good clusters from NB
    I = find(ptinx);
    cellids = CELLIDLIST(I);
else
    if isnumeric(I)  % I can be index set or list of cell IDs
        cellids = CELLIDLIST(I);
    else
        cellids = I;
        if ischar(cellids)
            cellids = {cellids};   % only one cell ID
        end
    end
end
cellids = cellids(:)';   % convert to row vector

% PSTH
[allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
    error_list] = poppsth(cellids,resdir,issave); %#ok<*ASGLU>
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
fnm = [resdir 'allPSTH.mat'];
if issave
    save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
        'tags','error_list')
end

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
list_name = ['error_list_' dsr];  % variable name
list_name = regexprep(list_name,'-','_');
list_name = regexprep(list_name,' ','_');
assignin('base',list_name,error_list)   % assign error list in base workspace
error_fname = fullfile(resdir,[list_name '.mat']);   % file name
save(error_fname,'error_list')   % save error list

% -------------------------------------------------------------------------
function [allpsth allpsth2 allspsth allspsth2 allstats ...
    ok_ids wn error_list] = poppsth(cellids,resdir,issave)

% Load CellBase
loadcb

% Time window
wn = [-2 2];   % in seconds
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
        alignevent = 'StimulusOn';
        
        % Calcualte PSTH
        [psth_shortRT, spsth_shortRT, spsth_se_shortRT, ~, spt_shortRT, stats_shortRT] = ...
            ultimate_psth(cellid,'trial',alignevent,wn,...
            'dt',dt,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','selectGoRT','filterinput',[0.15 0.5],'maxtrialno',Inf,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        [psth_longRT, spsth_longRT, spsth_se_longRT, ~, spt_longRT, stats_longRT] = ...
            ultimate_psth(cellid,'trial',alignevent,wn,...
            'dt',dt,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','selectGoRT','filterinput',[0.5 0.85],'maxtrialno',Inf,...
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
        
        % Save figure
        H = gcf;
        cellidt = regexprep(cellid,'\.','_');
        if issave
            fnm = [resdir cellidt '_PSTH.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot
%         Hr = rasterplot(spt);
%         if issave
%             fnm = [resdir cellidt '_RASTER.fig'];
%             saveas(Hr,fnm)
%         end
%         close(Hr)
        
        % Raster plot with sorted trials
        Hr2 = figure;
        viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','StimulusOff',...
            'eventtype','behav','ShowEvents',{{'StimulusOff'}},...
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
        fnm = [resdir cellidt '_SI.jpg'];   % save
%         saveas(Hr2,fnm)
        set(Hr2,'PaperPositionMode','auto')
        set(Hr2,'InvertHardcopy','off')
        print(Hr2,'-djpeg',fnm)
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

% -------------------------------------------------------------------------
function I = isnb(area1,area2)

nbareas = {'GP','GP/SI','SI','IC','RT/IC','EP','EA','EAC'};   % areas considered part of the basal forebrain (BF)
I = ismember(area1,nbareas);  % if 'primary' area is part of the BF

% -------------------------------------------------------------------------
function alignevent = findAlignEvent(cellid)

% Checking whether 'DeliverFeedback' event is available
sesstype = getvalue('session_type',cellid);
if isequal(sesstype,{'feedbackdelay'})
    alignevent = 'DeliverFeedback';
else
    alignevent = 'LeftWaterValveOn';
end