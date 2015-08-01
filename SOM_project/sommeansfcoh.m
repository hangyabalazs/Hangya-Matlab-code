function sommeansfcoh(cellids)
%SOMMEANSFCOH   Mean spike-field coherence.
%   SOMMEANSFCOH calculates average spike-field coherence between local
%   field potential (LFP) and single cell spiking data. The output files of
%   SOMSFCOH are used (see SOMSFCOH for further details).
%
%   SOMMEANSFCOH(I) uses cells defined by the index set I.
%   SOMMEANSFCOH(CELLIDS) uses cells defined by the list of CellIDs
%   (CELLIDS).
%
%   See also SOMSFCOH.

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
inpdir = [DATAPATH 'SOM\sfcoh_pv\'];

% Load spike-field coherence
NumCells = length(cellids);
allf = [];
allCxy = [];
for iC = 1:NumCells
    cellid = cellids{iC};   % current cell
    
    cellidt = regexprep(cellid,'\.','_');
    fn = [inpdir 'SFCOH_' cellidt '.mat'];
    load(fn)   % load coherence
    
    allf = [allf; fc']; %#ok<AGROW> % frequency vectors
    allCxy = [allCxy; Cxy']; %#ok<AGROW>  % coherence
end

% Check for compatibilty of frequency vectors
if any(any(allf-repmat(mean(allf),NumCells,1)))
    error('Frequency vector incompatibility.')
end

% Mean coherence
figure
plot(fc,mean(allCxy))
xlabel('Frequency (Hz)')
ylabel('Coherence')