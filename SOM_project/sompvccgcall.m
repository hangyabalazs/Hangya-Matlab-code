function sompvccgcall(cell_type)
%SOMPVCCGCALL   Calls SOMCCG_CONF_FILTER.
%   Caller for SOMCCG_CONF_FILTER. Copyright to Duda.
%
%   See also SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

% Handle input
% cell_type = pv_cells_v1;
CellVect = cell(1,length(cell_type));
CellPairs = cell(1,length(cell_type));
del_ind = zeros(1,length(CellPairs)); 

% Determine time window
sr = 1000;      % sampling rate
wn = 30 * sr / 1000;    % 2 * 30 ms window

% Clean cell pair list
for i = 1:length(CellPairs);
    CellVect{1,i} = find(strncmp(cell_type(i),cell_type,12)==1);
    CellPairs{1,i} = cell_type(CellVect{1,i});
    if length(CellPairs{1,i}) == 1
        del_ind(1,i) = i;
    end
end
del_ind = nonzeros(del_ind);
del_ind = del_ind';
CellPairs(del_ind) = [];
keep_cell = zeros(length(CellPairs),length(CellPairs));

for k = 1:length(CellPairs)-1
    for w = k+1:length(CellPairs)
        keep_cell(k,w) = strncmp(CellPairs{1,k}(1),CellPairs{1,w}(1),12);
    end
end
UniqueCellPairs = CellPairs(sum(keep_cell')==0);
counter = 1;
SpL = 50000;

% Preallocate
CCRMatrix = zeros(1000,2*wn+1);
LCCRMatrix = zeros(1000,2*wn+1);
UCCRMatrix = zeros(1000,2*wn+1);
MeanH0 = nan(1,1000);
SDH0 = nan(1,1000);
PairOfCells = cell(1,500);

% Calculate cross-correlation
for iCell = 1:length(UniqueCellPairs)
    Cells = combnk(UniqueCellPairs{1,iCell},2);
    [r c] = size(Cells);
    for iPair = 1:r
        Cell1 = Cells{iPair,1};
        Cell2 = Cells{iPair,2};
        tseg = findSegs2(Cell1);
        [ncc1,seltsind,selisi] = extractSegSpikes(Cell1,tseg);
        [ncc2,seltsind,selisi] = extractSegSpikes(Cell2,tseg);
        mn = min(ncc1(1),ncc2(1));  % only relative spike times count; avoid out of memory
        ncc1 = ncc1 - mn;
        ncc2 = ncc2 - mn;

        if length(ncc1) > SpL
            ncc1 = ncc1(1:SpL);   % include max 50000 spikes
        end
        if length(ncc2) > SpL
            ncc2 = ncc2(1:SpL);
        end
        if length(ncc1) > 100 && length(ncc2) > 100      % minimum 100 spikes
            PairOfCells{1,counter} = {Cell1 Cell2};
            clear Cell1 Cell2
            [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn);    % 1->2;
            CCRMatrix(counter,:) = ccr;
            LCCRMatrix(counter,:) = lwr;
            UCCRMatrix(counter,:) = upr;
            MeanH0(counter) = mean(mean(rccg,2),1);
            SDH0(counter) = mean(std(rccg,[],2),1);
            counter = counter + 1;
            disp('running......')
        end
    end
end

% Save
CCR = CCRMatrix(any(CCRMatrix,2),:);
LCCR = LCCRMatrix(any(LCCRMatrix,2),:);
UCCR = UCCRMatrix(any(UCCRMatrix,2),:);
MeanH0 = MeanH0(any(CCRMatrix,2));
SDH0 = SDH0(any(CCRMatrix,2));
% save pv_pv3 CCR UCCR LCCR PairOfCells CellPairs UniqueCellPairs MeanH0 SDH0
% save pv_pv_H0 MeanH0 SDH0
keyboard