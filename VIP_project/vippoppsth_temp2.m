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

