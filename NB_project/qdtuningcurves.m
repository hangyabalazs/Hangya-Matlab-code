function qdtuningcurves(I,issave)
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
resdir = [DATAPATH 'NB' fs 'tuningcurves2' fs];
 
% List of cellIDs
if isempty(I)
%     Lratio = getvalue('Lr_PC');
%     ID = getvalue('ID_PC');
%     area1 = getvalue('Area1');
%     area2 = getvalue('Area2');
%     inb = isnb(area1,area2);
%     ptinx = ID > 20 & Lratio < 0.15 & inb;
%     I = find(ptinx);
    selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
    ChAT = selectcell(selstr);
    I = cellfun(@findcellpos,ChAT);
end
I = I(:)';   % convert to row vector

% PSTH
poppsth(I,resdir,issave); %#ok<NASGU,*ASGLU>
% time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
% fnm = [resdir 'allPSTH.mat'];
% if issave
%     save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
%         'tags','problem_ids','problem_msg')
% end

% -------------------------------------------------------------------------
function poppsth(I,resdir,issave)

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
        [psth, spsth, ~, tags, spt] = ultimate_psth(cellid,'trial',...
            'LeftWaterValveOn',wn,...
            'dt',0.001,'display',true,'sigma',0.02,'parts','#StimulusDuration',...
            'isadaptive',0,'event_filter','custom','filterinput','Hit==1',...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        
        % Concatenate data from different cells
%         allpsth = [allpsth; psth]; %#ok<AGROW>
%         allspsth = [allspsth; spsth]; %#ok<AGROW>
%         allstats = [allstats; stats]; %#ok<AGROW>
%         ok_ids{end+1} = cellid; %#ok<AGROW>
        
        % Plot
        figure
        plot(spsth')
        legend(tags)

        % Save figure
        if issave
            H = gcf;
            cellidt = regexprep(cellid,'\.','_');
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
        
        
    catch ME
        
        % Error handling
        problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
        problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
    end
end

% Standardize
% allpsth2 = allpsth;
% allspsth2 = allspsth;
% for k = 1:size(allpsth,1)
%     allpsth2(k,:) = standardize(allpsth(k,:));
%     allspsth2(k,:) = standardize(allspsth(k,:));
% end

% -------------------------------------------------------------------------
function I = isnb(area1,area2)

nbareas = {'GP','GP/SI','SI','IC','RT/IC','EP','EA','EAC'};   % aeas considered part of the basal forebrain (BF)
I = ismember(area1,nbareas);  % if 'primary' area is part of the BF