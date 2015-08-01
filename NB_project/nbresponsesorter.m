function error_list = nbresponsesorter(I,issave)
%NBRESPONSESORTER   Peri-event time histogram.
%   NBRESPONSESORTER(I,ISSAVE) calculates 'doubly adaptive' PSTHs for a
%   set of cells (see DAPSTH) aligned to 'Go' or 'No-go' response onset
%   ('Hit' or 'False Alarm'). Statistical tests are performed to probe
%   significant firing rate changes after responses (see PSTH_STATS). Input
%   parameters:
%       I - list of cell IDs or index set to CELLIDLIST  (see CellBase 
%           documentation); if empty or not specified, all well-separated
%           cells are selected (ID>20, L-ratio<0.15; see LRATIO) from basal
%           forebrain areas
%       ISSAVE - controls saving
%
%   ERROR_LIST = NBRESPONSESORTER(I,ISSAVE) returns a structure with all
%   caught errors. ERROR_LIST includes cell IDs, error messages and the
%   captured exception variables. A time stamped ERROR_LIST is also
%   assigned in base workspace and saved to the results directory
%   automatically.
%
%   See also DAPSTH, ULTIMATE_PSTH, PSTH_STATS, LRATIO, GAP_STATISTICS and
%   NBPOPPSTH_CALL.

%   Edit log: BH 7/3/12, 8/13/12, 4/23/12
 
% Pass the control to the user in case of error
dbstop if error
choosecb('NB')

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
resdir = [DATAPATH 'NB' fs 'responsesorter_newdata2' fs];   % results directory

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
function [allpsth allpsth2 allspsth allspsth2 allstats ok_ids wn...
    error_list] = poppsth(cellids,resdir,issave)

% Load CellBase
loadcb

% Time window
wn = [-0.6 0.6];   % in seconds

% Call 'main'
% Preallocate
allpsth = [];
allspsth = [];
allstats = struct([]);
error_list = struct('cellid',{},'message',{},'exception',{});  % keep a list of caught errors 
errinx = 0;
ok_ids = {};
NumCells = length(cellids);
for k = 1:NumCells
    cellid = cellids{k};   % cell ID
    disp(cellid)
    
    if strcmp(getvalue('session_type',cellid),'pavlovian')
        disp([ cellid ' was recorded in Pavlovian rask.'])
        continue
    end
    try
        
        % Check for 'DeliverFeedback' event
        alignfilter = 'FalseAlarm==1';
        alignevent = findAlignEvent(cellid);
        
        % Calcualte PSTH
        [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
            alignevent,wn,...
            'dt',0.001,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        
        % Concatenate data from different cells
        allpsth = [allpsth; psth]; %#ok<AGROW>
        allspsth = [allspsth; spsth]; %#ok<AGROW>
        allstats = [allstats; stats]; %#ok<AGROW>
        ok_ids{end+1} = cellid; %#ok<AGROW>
        
        % Save figure
%         if issave
%             H = gcf;
%             cellidt = regexprep(cellid,'\.','_');
%             fnm = [resdir cellidt '_PSTH.fig'];
%             saveas(H,fnm)
%         end
%         close(H)
%         
%         % Raster plot
%         Hr = rasterplot(spt);
%         if issave
%             fnm = [resdir cellidt '_RASTER.fig'];
%             saveas(Hr,fnm)
%         end
%         close(Hr)
        
        
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
    alignevent = 'LeftPortIn';
end