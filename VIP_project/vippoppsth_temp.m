function vippoppsth

% load

load('C:\Balazs\_analysis\VIP\vipisinfluenced11\validity.mat')
load('C:\Balazs\_analysis\VIP\vipisinfluenced11\p_val.mat')

% with FR criterion

inx_act = logical(vldty') & (p_act<0.015) & (baseline>2);
inx_act(105) = 1;   % tagged
inx_inh = logical(vldty') & (p_inh<0.015) & (baseline>2);

tagged = find(inx_act&activation_start<3);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = ai(inx);
inhibited_activated = ai(~inx);
activated = [activated activated_inhibited];
inhibited = [inhibited inhibited_activated];

csb = whichcb;
if ~isequal(csb,'VIP')
    disp('Switch to VIP CellBase!')
    return
end
loadcb

% Plot WS ramp neurons poppsth
BurstPSTH='off';
Sort= 'none';
Normalization = 'zscore'; %'zscore'; %zscore, none etc.
TriggerName='BurstOn';  %'NextForwardRun3';
Window = [-0.02 0.06];
FigNum = 40;
dt = 0.001;
sigma  = 0.05;
Overlay = 'on';
PlotDashedEvent = '';
PlotDashedCondition = '';


handles3=viewpoppsth3_stim(CELLIDLIST(activated),[],'TriggerName',TriggerName,'Normalization',Normalization,...
'window',Window,'FigureNum',FigNum,'Overlay',Overlay,'BurstPSTH',BurstPSTH,...
'dt',dt,'sigma',sigma,...
'PlotDashedEvent',PlotDashedEvent,'PlotDashedCondition',PlotDashedCondition);


viewpoppsth3_stim(CELLIDLIST(activated))