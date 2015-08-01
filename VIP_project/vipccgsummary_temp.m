%%

% Load cross-correlograms
load('C:\Balazs\_analysis\VIP\CCG_A1\CCG_matrices.mat')


%%

% Load PSTH variables
load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')

% Groups
frlim = 1;   % lower firing rate limit for detecting inhibition
tagged = logical(vldty) & (isact==2);   % tagged cells
inx_act = logical(vldty) & (isact==1) & (isact~=2);  % activated cells
inx_inh = logical(vldty) & (isinh) & (baseline>frlim) & (isact~=2);   % inhibited cells; firing rate > 1Hz criterion
activated = find(inx_act&~inx_inh);  % indices of activated only cells
inhibited = find(inx_inh&~inx_act);  % indices of inhibited only cells
ai = find(inx_act&inx_inh);   % indices of cells with dual effect
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));  % activated, then inhibited
inhibited_activated = ai(~inx);   % inhibited, then activated
activated = [activated; activated_inhibited];   % activated and activated-inhibited
inhibited = [inhibited; inhibited_activated];   % inhibited and inhibited-activated
tagged = find(tagged);   % indices of tagged cells
[~, iainx] = intersect(inhibited,inhibited_activated);   % indices of inh-act within inh

%% 

% Load CellBase
loadcb   % get list of cell IDs

% Find pairs with inhibited cells
inhinx = cellfun(@(s)ismember(s,CELLIDLIST(inhibited)),PairOfCells);   % indicats whether a cell ID in PairOfCells corresponds to an inhibited neuron
inhibited_poc = any(inhinx,2);   % pairs including inhibited cells
inh_poc_inx = find(inhibited_poc);   % convert to indices to CCR matrices

% Plot inhibited CCGs with significant effects
LCCR2 = LCCR;
LCCR2(LCCR2<1) = -0.1;   % don't consider fluctuations between 0 and 1 (non-stationarity)
NumInh = length(inh_poc_inx);
mins = nan(NumInh,1);  % effect size for inhibition
maxs = nan(NumInh,1);  % effect size for excitation
for k = 1:NumInh   % loop through CCGs with inhibited cells 
    k2 = inh_poc_inx(k);   % current index to CCG matrices
    if isequal(inhinx(k2,:),[0 1])   % flip so that the inhibited cell is first
        ccr = fliplr(CCR(k2,:));
        lccr = fliplr(LCCR2(k2,:));
        uccr = fliplr(UCCR(k2,:));
    else
        ccr = CCR(k2,:);
        lccr = LCCR2(k2,:);
        uccr = UCCR(k2,:);
    end
    if any(ccr(27:35)<ccr(27:35)) || ...
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

%%

% Find pairs with activated cells
actinx = cellfun(@(s)ismember(s,CELLIDLIST(activated)),PairOfCells);   % indicats whether a cell ID in PairOfCells corresponds to an activated neuron
activated_poc = any(actinx,2);   % pairs including activated cells
act_poc_inx = find(activated_poc);   % convert to indices to CCR matrices

% Plot activated CCGs with significant effects
LCCR2 = LCCR;
LCCR2(LCCR2<1) = -0.1;   % don't consider fluctuations between 0 and 1 (non-stationarity)
NumAct = length(act_poc_inx);
mins = nan(NumAct,1);  % effect size for inhibition
maxs = nan(NumAct,1);  % effect size for excitation
for k = 1:NumAct   % loop through CCGs with activated cells 
    k2 = act_poc_inx(k);   % current index to CCG matrices
    if isequal(actinx(k2,:),[0 1])   % flip so that the activated cell is first
        ccr = fliplr(CCR(k2,:));
        lccr = fliplr(LCCR2(k2,:));
        uccr = fliplr(UCCR(k2,:));
    else
        ccr = CCR(k2,:);
        lccr = LCCR2(k2,:);
        uccr = UCCR(k2,:);
    end
    if any(ccr(27:35)<ccr(27:35)) || ...
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