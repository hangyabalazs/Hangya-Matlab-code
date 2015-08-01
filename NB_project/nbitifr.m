function nbitifr(cellids,issave)
%NBITIFR   Firing during the foreperiod.
%   NBITIFR(CELLIDS,ISSAVE) plots raster plots and PSTHs aligned to trial
%   starts to help the evaluation of potential firing rate changes based on
%   temporal expectations. Spikes are prealigned to the ITI event (from
%   LastITIBegins to LastITIEnds). Two versions are saved: partioned
%   according to outcome or not partitioned.
%   Input parameters: 
%       CELLIDS - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, cholinergic cells
%           are selected
%       ISSAVE - controls saving
%
%   See also VIEWCELL2B.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   30-Oct-2013

%   Edit log: BH 10/30/13

% Pass the control to the user in case of error
dbstop if error

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1 || isempty(cellids)
%     loadcb   % load CellBase
%     cellids = CELLIDLIST;
    
    % All ChAT and pChAT cells
    ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
    pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
    allChAT = [ChAT pChAT];
    cellids = allChAT;
else
    if isnumeric(cellids)
        loadcb   % load CellBase
        cellids = CELLIDLIST(cellids);   % index set to CELLIDLIST
    elseif ischar(cellids)
        cellids = {cellids};   % only one cellID passed
    elseif iscellstr(cellids)
        % list of cell IDs
    else
        error('nbitifr:inputArg','Unsupported format for cell IDs.')
    end
end
NumCells = length(cellids);   % number of cells

% Directories
global DATAPATH
% resdir = [DATAPATH 'NB\itifr\'];
resdir = [DATAPATH 'NB\attentioncells\'];

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running NBITIFR...');  % progress indicator
global WB
WB(end+1) = wb;

% PSTH aligned to trial starts
for iC = 1:NumCells   % loop through cells
    cellid = cellids{iC};   % current cell
    TE = loadcb(cellid,'TrialEvents');  % load trial events
    
    % Add new event for ITI
    TS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
    if isequal(findcellstr(TS.events(:,1),'ITI'),0)   % prealign spikes to trial start if it has not happened before
        prealignSpikes(cellid,'FUNdefineEventsEpochs',...
            @defineEventsEpochs_iti,'filetype','behav',...
            'ifsave',1,'ifappend',1)
    end
    
    % Raster and PSTH aligned to trial start
    H1 = figure;
    viewcell2b(cellid,'TriggerName','ITI','SortEvent','TrialStart','eventtype','behav',...
        'LastEvents','LastITIEnds','ShowEvents',{{'StimulusOn'}},'Partitions','#ResponseType','window',[-2 5])
    maximize_figure(H1)   % partitioned on outcome
    
    H2 = figure;
    viewcell2b(cellid,'TriggerName','ITI','SortEvent','TrialStart','eventtype','behav',...
        'LastEvents','LastITIEnds','ShowEvents',{{'StimulusOn'}},'Partitions','all','window',[-2 5])
    maximize_figure(H2)   % not partitioned
    
    % Save
    if issave
        cellidt = regexprep(cellid,'\.','_');
        fnm = [resdir cellidt '_ITI1.jpg'];
        saveas(H1,fnm)
        fnm = [resdir cellidt '_ITI2.jpg'];
        saveas(H2,fnm)
    end
    
    close all
    waitbar(iC/NumCells)
end
close(wb)