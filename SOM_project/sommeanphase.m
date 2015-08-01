function sommeanphase(cellids)
%SOMMEANPHASE   Summary phase statistics.
%   SOMMEANPHASE plots mean phase and mean vectors in polar coordinates.
%   The output files of SOMPHASEHIST are used (see SOMPHASEHIST for further
%   details).
%
%   SOMMEANPHASE(I) uses cells defined by the index set I.
%   SOMMEANPHASE(CELLIDS) uses cells defined by the list of CellIDs
%   (CELLIDS).
%
%   See also SOMPHASEHIST.

% Input argument check
if nargin < 1
    loadcb   % load CellBase
    cellids = CELLIDLIST; 
else
    if isnumeric(cellids)
        loadcb   % load CellBase
        cellids = CELLIDLIST(cellids);
    end
end
if ischar(cellids)   % convert 'cellids' to cell
    cellids = {cellids};
end

% Directories
global DATAPATH
inpdir = [DATAPATH 'SOM\phasehist_som_theta\'];

% Load phase statistics
NumCells = length(cellids);
allang = [];
allmvl = [];
figure
for iC = 1:NumCells
    cellid = cellids{iC};   % current cell
    
    cellidt = regexprep(cellid,'\.','_');
    fn = [inpdir 'PHASE_' cellidt '.mat'];
    load(fn)   % load phase
    
    allang = [allang hang]; %#ok<AGROW> % mean angle
    allmvl = [allmvl hmvl]; %#ok<AGROW>  % mean vector length
    
    % Plot
    CMP = compass(real(0.15*exp(1).^(1i*0.15)),imag(0.15*exp(1).^(1i*0.15)));
    set(CMP,'Color','white')
    hold on
    CMP = compass(real(ftm),imag(ftm));
    set(CMP,'Color','blue','LineWidth',2)
end