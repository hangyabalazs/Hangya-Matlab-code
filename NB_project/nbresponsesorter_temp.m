function nbresponsesorter(I,issave)
%NBRESPONSESORTER   Cluster analysis on PSTHs.
%   NBRESPONSESORTER(I,ISSAVE) calculates 'doubly adaptive' PSTHs for a set
%   of cells (see DAPSTH) aligned to 'Go' response onset ('Hit').
%   Statistical tests are performed to probe significant firing rate
%   changes after responses (see PSTH_STATS). The PSTHs are then clustered
%   into response categories and natural clustering is tested using
%   GAP_STATISTICS. Input parameters: 
%       I - index set to CELLIDLIST (see CellBase documentation); if empty
%           or not specified, all well-separated cells are selected (ID>20,
%           L-ratio<0.15; see LRATIO) from basal forebrain areas.
%       ISSAVE - controls saving
%
%   See also DAPSTH, ULTIMATE_PSTH, PSTH_STATS, LRATIO, GAP_STATISTICS and
%   NBPOPPSTH_CALL.

%   Edit log: BH 7/3/12, 8/13/12
 
% Pass the control to the user in case of error
dbstop if error
choosecb('NB')

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
resdir = [DATAPATH 'NB' fs 'responsesorter6' fs];
 
% List of cellIDs
if isempty(I)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    vldty = getvalue('validity');
    area1 = getvalue('Area1');
    area2 = getvalue('Area2');
    inb = isnb(area1,area2);
    ptinx = vldty == 1 & ID > 20 & Lratio < 0.15 & inb;
    I = find(ptinx);
end
I = I(:)';   % convert to row vector

% PSTH
[allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
    problem_ids problem_msg] = poppsth(I,resdir,issave); %#ok<NASGU,*ASGLU>
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
fnm = [resdir 'allPSTH.mat'];
if issave
    save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
        'tags','problem_ids','problem_msg')
end

% Clustering responses
dec = 50;
naninx = all(isnan(allpsth),2);
allpsth(naninx,:) = [];
[khat Gap s_k C] = gap_statistics(allpsth,dec);
fnm = [resdir 'gapstat.mat'];
if issave
    save(fnm,'khat','Gap','s_k','C')
end

% Visualization
% figure
% for k = 1:khat
%     eval(['subplot(' num2str(khat) ',1,' num2str(k) ')'])
%     imagesc(allspsth(C==k,:))
%     set(gca,'XTick',[])
%     colorbar
% end

% -------------------------------------------------------------------------
function [allpsth allpsth2 allspsth allspsth2 allstats ok_ids wn...
    problem_ids problem_msg] = poppsth(I,resdir,issave)

% Load CellBase
loadcb

% Time window
wn = [-0.6 0.6];   % in seconds

% Call 'main'
% Preallocate
allpsth = [];
allspsth = [];
allstats = struct([]);
problem_ids = {};
problem_msg = {};
ok_ids = {};
for k = I
    cellid = CELLIDLIST{k}; %#ok<USENS>
    disp(cellid)
    try
        
        % Calcualte PSTH
        [psth, spsth, ~, spt, stats] = ultimate_psth(cellid,'trial',...
            'LeftWaterValveOn',wn,...
            'dt',0.001,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','custom','filterinput','Hit==1','maxtrialno',Inf,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        
        % Concatenate data from different cells
        allpsth = [allpsth; psth]; %#ok<AGROW>
        allspsth = [allspsth; spsth]; %#ok<AGROW>
        allstats = [allstats; stats]; %#ok<AGROW>
        ok_ids{end+1} = cellid; %#ok<AGROW>
        
        % Save figure
        if issave
            H = gcf;
            cellidt = regexprep(cellid,'\.','_');
            fnm = [resdir cellidt '_PSTH.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot
        Hr = rasterplot(spt);
        if issave
            fnm = [resdir cellidt '_RASTER.fig'];
            saveas(Hr,fnm)
        end
        close(Hr)
        
        
    catch ME
        
        % Error handling
        problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
        problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
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