function error_list = vipresponsesorter(I,issave)
%VIPBRESPONSESORTER   Peri-event time histogram.
%   VIPRESPONSESORTER(I,ISSAVE) calculates 'doubly adaptive' PSTHs for a
%   set of cells (see DAPSTH) aligned to 'Go' or 'No-go' response onset
%   ('Hit' or 'False Alarm'). Statistical tests are performed to probe
%   significant firing rate changes after responses (see PSTH_STATS). Input
%   parameters: 
%       I - list of cell IDs or index set to CELLIDLIST  (see CellBase 
%           documentation); if empty or not specified, all cells are
%           selected
%       ISSAVE - controls saving
%
%   ERROR_LIST = VIPRESPONSESORTER(I,ISSAVE) returns a structure with all
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
choosecb('VIP_gonogo')
trialtype = 'CR';
align = 'tone';

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1
    I = [];
end
 
% Directories
global DATAPATH
fs = filesep;
if isequal(whichcb,'VIP_gonogo')
    resdir = [DATAPATH 'VIP' fs 'responsesorter2_' trialtype '_newwin4_stimaligned' fs];   % results directory
elseif isequal(whichcb,'VIP_gonogo2')
    resdir = [DATAPATH 'VIP' fs 'responsesorter2_' trialtype '_newwin4_stimaligned_shock' fs];   % results directory
else
    error('vipresponsesorter:cellBase','Unsupported CellBase.')
end
        
% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% List of cellIDs
if isempty(I)
    cellids = CELLIDLIST;
else
    if isnumeric(I)  % I can be index set or list of cell IDs
        cellids = CELLIDLIST(I);
    elseif ischar(I)  % one cell ID is passed
        cellids = {I};
    else
        cellids = I;
    end
end
cellids = cellids(:)';   % convert to row vector

% PSTH
[allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
    error_list] = poppsth(cellids,resdir,issave,trialtype,align); %#ok<*ASGLU>
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
    error_list] = poppsth(cellids,resdir,issave,trialtype,align)

% Load CellBase
loadcb

% Time window
% wn = [-0.6 0.6];   % in seconds
wn = [-1.2 2.2];   % in seconds

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
    
    try
        
        % Align to 'Lick' event
        switch trialtype
            case 'Hit'
                alignfilter = 'ResponseType==1';
            case 'FA'
                alignfilter = 'ResponseType==2';
            case 'CR'
                alignfilter = 'ResponseType==3';
            case 'Miss'
                alignfilter = 'ResponseType==4';
            otherwise
                error('vipresponsesorter:trialType','Invalid trial type.')
        end
        if isequal(whichcb,'VIP_gonogo')
            switch align
                case 'tone'
                    alignevent = 'StimOnset';   % VIP_gonogo CellBase
                case 'fb'
                    alignevent = 'Lick';   % VIP_gonogo CellBase
                otherwise
                    error('vipresponsesorter:alignType','Invalid align event.')
            end
        elseif isequal(whichcb,'VIP_gonogo2')
            switch align
                case 'tone'
                    alignevent = 'StimulusOn';   % VIP_gonogo2 CellBase
                case 'fb'
                    switch align
                        case 'Hit'
                            alignevent = 'LeftWaterValveOn';   % for Hit, VIP_gonogo2 CellBase
                        case 'FA'
                            alignevent = 'LeftPortIn';   % for FA, VIP_gonogo2 CellBase
                    end
                otherwise
                    error('vipresponsesorter:alignType','Invalid align event.')
            end
        else
            error('vipresponsesorter:cellBase','Unsupported CellBase.')
        end
        
        % Calcualte PSTH
        [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
            alignevent,wn,...
            'dt',0.001,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
            'baselinewin',[-1 -0.6],'testwin',[0 0.4],'relative_threshold',0.1);
        
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
        % Concatenate data from different cells
        allpsth = [allpsth; nan(1,size(allpsth,2))]; %#ok<AGROW>
        allspsth = [allspsth; nan(1,size(allpsth,2))]; %#ok<AGROW>
        fld = fieldnames(allstats);   % fieldnames
        str = [fld'; repmat({''',NaN,'''},size(fld'))];
        str = str(:)';
        str = cell2mat(str);
        str = ['''' str(1:end-2)];
        str = ['struct(' str ');'];
        allstats = [allstats; eval(str)]; %#ok<AGROW>
        
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