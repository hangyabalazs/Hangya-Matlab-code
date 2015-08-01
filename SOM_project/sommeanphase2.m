function sommeanphase2(cellids)
%SOMMEANPHASE2   Summary phase statistics.
%   SOMMEANPHASE2 plots mean phase and mean vectors in polar coordinates.
%   Only significantly phase locked cells are included (p<0.05,
%   Rayleigh-test) The output files of SOMPHASEHIST are used (see
%   SOMPHASEHIST for further details).
%
%   SOMMEANPHASE2(I) uses cells defined by the index set I.
%   SOMMEANPHASE2(CELLIDS) uses cells defined by the list of CellIDs
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
allang_sig = [];
allmvl = [];
figure
cntr = 0;
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
    if p_rayleigh < 0.05
        allang_sig = [allang_sig hang]; %#ok<AGROW>
        CMP = compass(real(ftm),imag(ftm));
        set(CMP,'Color','red','LineWidth',2)
        cntr = cntr + 1;
    end
end

% Histogram of preferred phase angles
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
nm = histc(rad2deg(allang_sig),edges);   % phase histogram
nm = nm(1:end-1);
figure
B = bar(cnts,nm);
set(B,'FaceColor','k')

keyboard