function vipccgsummary

% Load cross-correlograms
load('C:\Balazs\_analysis\VIP\CCG_mPFC\CCG_matrices.mat')
LCCR2 = LCCR;
LCCR2(LCCR2<1) = -0.1;   % don't consider fluctuations between 0 and 1 (non-stationarity)

% Load PSTH variables
load('C:\Balazs\_analysis\VIP\psth_variables.mat')

% Load CellBase
loadcb   % get list of cell IDs

% Groups
tagged =  logical(vldty) & (isact==2);
inx_act = logical(vldty) & (isact==1) & (baseline>2);
inx_inh = logical(vldty) & (isinh) & (baseline>2) & (isact~=2);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort([ai(inx); 10]);    % 10 is activated-inhibited
inhibited_activated = setdiff(ai(~inx),10);
activated = [activated; activated_inhibited];
inhibited = [inhibited; inhibited_activated];
inhibited_short = intersect(inhibited,find(inhibition_peak<40));
inhibited_long = intersect(inhibited,find(inhibition_peak>=40));
inhibited = inhibited_short;

% Find pairs with inhibited cells
inhinx = cellfun(@(s)ismember(s,CELLIDLIST(inhibited)),PairOfCells);   % indicates whether a cell ID in PairOfCells corresponds to an inhibited neuron
inhibited_poc = any(inhinx,2);   % pairs including inhibited cells
inh_poc_inx = find(inhibited_poc);   % convert to indices to CCR matrices

% Plot CCGs for inhibted cells
plotccgs(CCR,LCCR2,UCCR,inh_poc_inx,inhinx)

% Find pairs with activated cells
actinx = cellfun(@(s)ismember(s,CELLIDLIST(activated)),PairOfCells);   % indicates whether a cell ID in PairOfCells corresponds to an activated neuron
activated_poc = any(actinx,2);   % pairs including activated cells
act_poc_inx = find(activated_poc);   % convert to indices to CCR matrices

% Plot CCGs for activated cells
plotccgs(CCR,LCCR2,UCCR,act_poc_inx,actinx)

% Find pairs with activated cells
taginx = cellfun(@(s)ismember(s,CELLIDLIST(tagged)),PairOfCells);   % indicates whether a cell ID in PairOfCells corresponds to a tagged neuron
tagged_poc = any(taginx,2);   % pairs including tagged cells
tag_poc_inx = find(tagged_poc);   % convert to indices to CCR matrices

% Plot CCGs for tagged cells
plotccgs(CCR,LCCR2,UCCR,tag_poc_inx,taginx)

keyboard

% -------------------------------------------------------------------------
function plotccgs(CCR,LCCR,UCCR,poc_inx,pinx)


% Plot CCGs with significant effects
NumIPairs = length(poc_inx);   % number of pairs having 'identified' cells
mins = nan(NumIPairs,1);  % effect size for inhibition
maxs = nan(NumIPairs,1);  % effect size for excitation
for k = 1:NumIPairs   % loop through CCGs for 'identified' cells
    k2 = poc_inx(k);   % current index to CCG matrices
    if isequal(pinx(k2,:),[0 1])   % flip so that the 'identified' cell is first
        ccr = fliplr(CCR(k2,:));
        lccr = fliplr(LCCR(k2,:));
        uccr = fliplr(UCCR(k2,:));
    else
        ccr = CCR(k2,:);
        lccr = LCCR(k2,:);
        uccr = UCCR(k2,:);
    end
    if any(ccr(27:35)<lccr(27:35)) || ...
            any(ccr(27:35)>uccr(27:35))   % effect within 4 ms from center (0)
        
        % Plot
        figure
        plot(ccr,'k')   % plot CCG
        hold on
        plot(lccr,'r')  % plot lower confidence interval
        plot(uccr,'r')  % plot upper confidence interval
        
        % Effect size
        minpc = - mean(ccr([1:26 36:61])) + min(ccr(27:35));   % diff. between mean and trough
        maxpc = max(ccr(27:35)) - mean(ccr([1:26 36:61]));   % diff. between peak and mean
        if any(ccr(27:35)<lccr(27:35))   % inhibition in CCG
            mins(k2) = minpc;  % effect size for inhibition
        end
        if any(ccr(27:35)>uccr(27:35))   % excitaion in CCG
            maxs(k2) = maxpc;  % effect size for excitation
        end
    end
end