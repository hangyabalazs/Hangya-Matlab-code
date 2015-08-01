function sompvccgcall_hetero(cell_type,test_cells)
%SOMPVCCGCALL_HETERO   Calls SOMCCG_CONF_FILTER for heterolog pairs.
%   Caller for SOMCCG_CONF_FILTER. Copyright to Duda.
%
%   See also SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

% Handle input
% cell_type = pv_cells_v1;
% test_cells = non_tagged;

% Determine time window
sr = 1000;      % sampling rate
wn = 30 * sr / 1000;    % 2 * 30 ms window

% Cell pair loop for CCG
method = 'random';
counter = 1;
limit_spikes = 50000;   % include max 50000 spikes
CCRMatrix = zeros(1000,61);
PairOfCells = cell(1,3000);
LCCRMatrix = zeros(1000,61);
UCCRMatrix = zeros(1000,61);
MeanH0 = nan(1,1000);
SDH0 = nan(1,1000);
BefEV = 2;
AftEV = 2;
for iPair = 1:length(cell_type)
    [r,s,t,u] = cellid2tags(cell_type(1,iPair));   % get cell pairs
    sess_cells = findcell('rat',r,'session',s);
    sess_cells = sess_cells(find(strncmp(r,sess_cells,4)==1));
    sess_cells = setdiff(sess_cells,cell_type);
    sess_cells = intersect(test_cells,sess_cells);
    for iCouple = 1:length(sess_cells)
        Cell1 = cell_type(1,iPair);
        Cell2 = sess_cells(iCouple);
        switch method
            case 'event_triggered'
                ST1 = loadcb(Cell1 ,'EVENTSPIKES');
                ST2 = loadcb(Cell2 ,'EVENTSPIKES');
                TE = loadcb(Cell1,'TrialEvents');
                BEvent = 'HomeZoneOut1';
                ind = find(strcmp(BEvent,ST1.events(:,1)));
                try ncc1 = ST1.event_stimes{ind}(1,:); catch; end
                try ncc2 = ST2.event_stimes{ind}(1,:); catch; end
                for iEV = 1:length(ncc1)
                   ncc1{1,iEV} = ncc1{1,iEV}(ncc1{1,iEV}>-BefEV & ncc1{1,iEV}<AftEV);
                   ncc2{1,iEV} = ncc2{1,iEV}(ncc2{1,iEV}>-BefEV & ncc2{1,iEV}<AftEV);
                   ncc1{1,iEV} = ncc1{1,iEV} + TE.HomeZoneOut1(iEV);
                   ncc2{1,iEV} = ncc2{1,iEV} + TE.HomeZoneOut1(iEV);
                end
                ncc1 = cell2mat(ncc1');
                ncc2 = cell2mat(ncc2');
            case 'random'
                try 
                    tseg = findSegs2(Cell1);
                    [ncc1,seltsind,selisi] = extractSegSpikes(Cell1,tseg);
                    [ncc2,seltsind,selisi] = extractSegSpikes(Cell2,tseg);
                catch
                    disp('can not extract the segment');
                end
        end
        if length(ncc1) > limit_spikes;
            ncc1 = ncc1(1:limit_spikes);
        end
        if length(ncc2) > limit_spikes
            ncc2=ncc2(1:limit_spikes);
        end
        if length(ncc1) > 100 && length(ncc2) > 100     % minimum 100 spikes
            PairOfCells{1,counter} = {Cell1 Cell2};
            clear Cell1 Cell2
            [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn);    % 1->2
            close all
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
% save pv_nontagged2 CCR UCCR LCCR PairOfCells MeanH0 SDH0
save pv_nontagged_H0 MeanH0 SDH0
keyboard