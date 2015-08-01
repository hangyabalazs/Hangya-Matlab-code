function nbplotpsth

% Pass the control to the user in case of error
dbstop if error
 
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
resdir = [DATAPATH 'NB' fs 'psth_summary' fs];
 
% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Find putative tagged cells
if isempty(I)
    mode = 'groups';
    
    % PSTH statistics
    allstats = getvalue('FA_psth_stats');
    allstats = nan2cellstruct(allstats);
    activation_peak = [allstats.activation_peak];   % peak time of activation
    Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
    inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
    Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
    
    % Groups of activated and inhibited cells
    inx_act = Wpa < 0.01;   % significant activation
    inx_inh = Wpi < 0.01;   % significant inhibition
    activated = find(inx_act&~inx_inh);    % indices for activated cells
    inhibited = find(inx_inh&~inx_act);    % indices for inhibited cells
    ai = find(inx_act&inx_inh);
    inx = activation_peak(ai) < inhibition_peak(ai);
    activated_inhibited = sort(ai(inx));
    inhibited_activated = ai(~inx);
    activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
    inhibited = [inhibited'; inhibited_activated'];
    
    % Identified NB cells
    selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
    ChAT = selectcell(selstr);   % cell IDs for ChAT cells
    selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
    pChAT = selectcell(selstr);  % cell IDs for putative ChAT cells
    allChAT = [ChAT pChAT];   % cell IDs for identified and putative ChAT cells
    selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
    NT = selectcell(selstr);   % cell IDs for ChAT cells
    
    % Final groups
    NTact = intersect(NT,CELLIDLIST(activated));   % activated NT cells
    NTinh = intersect(NT,CELLIDLIST(inhibited));   % inhibited NT cells
    
else
    mode = 'custom';
    if isnumeric(I)
        cellids = CELLIDLIST(I);   % index set to CELLIDLIST
    elseif ischar(I)
        cellids = {I};   % only one cellID passed
    elseif iscellstr(I)
        cellids = I;   % list of cell IDs
    else
        error('nbplotpsth:inputArg','Unsupported format for cell IDs.')
    end
end

% PSTH
switch mode
    case 'group'   % identified and unindentified groups
        main(allChAT)
        main(NTact)
        main(NTinh)
    case 'custom'   % custom-list of cell IDs
        main(cellids)
end