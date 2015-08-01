%%

load('EVENTS.mat')

%%

load('F:\NB_cellbase\n027\120215e\2012-02-15_21-07-17_LaserStimProtocol2.mat')

%%

inx = find(Events_Nttls==16384);
stim_nttl_ts = Events_TimeStamps(inx);
burststart_inx = [0 find(diff(stim_nttl_ts)>2)];

%%

lb =length(BurstDur);
burststart_inx(end+1) = length(inx);
prec = 0.001;
Events_Nttls_new = [Events_Nttls(1:inx(1)-1) 1];
Events_TimeStamps_new = [Events_TimeStamps(1:inx(1)-1) Events_TimeStamps(inx(1))-2*prec];
Events_EventStrings_new = [Events_EventStrings(1:inx(1)-1); {'C:\Documents and Settings\Administrator\Desktop\nb027\2012-02-15_21-07-17_LaserStimProtocol2'}];
Events_EventIDs_new = [Events_EventIDs(1:inx(1)-1) 3];
for k = 1:lb
    Events_Nttls_new = ...
        [Events_Nttls_new 4 Events_Nttls((inx(1)+burststart_inx(k)*2):...
        (inx(1)+burststart_inx(k+1)*2)-1)];
    str = {['NTrial=' num2str(k) ';' ...
        'BurstPulseDur=' num2str(BurstPulseDur(k)) ';' ...
        'BurstPulseIPI=' num2str(BurstPulseIPI(k)) ';' ...
        'BurstNPulse=' num2str(BurstNPulse(k)) ';' ...
        'BurstPulsePower=' num2str(BurstPulsePower(k)) ';']};
    Events_EventStrings_new = [Events_EventStrings_new; str; ...
        Events_EventStrings((inx(1)+burststart_inx(k)*2):...
        (inx(1)+burststart_inx(k+1)*2)-1)];
    Events_TimeStamps_new = [Events_TimeStamps_new Events_TimeStamps(inx(1)+burststart_inx(k)*2)-prec Events_TimeStamps((inx(1)+burststart_inx(k)*2):...
        (inx(1)+burststart_inx(k+1)*2)-1)];
    Events_EventIDs_new = [Events_EventIDs_new 3 ...
        Events_EventIDs((inx(1)+burststart_inx(k)*2):...
        (inx(1)+burststart_inx(k+1)*2)-1)];
end


%% 

Events_TimeStamps = Events_TimeStamps_new;
Events_Nttls = Events_Nttls_new;
Events_EventStrings = Events_EventStrings_new;
Events_EventIDs = Events_EventIDs_new;