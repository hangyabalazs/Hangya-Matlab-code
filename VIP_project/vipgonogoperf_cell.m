%% Performance

lastinx = max(find(Hit==1|FalseAlarm==1));

GoPerf = sum(Hit(1:lastinx)==1) / (sum(Hit(1:lastinx)==1|Miss(1:lastinx)==1));
NoGoPerf = sum(CorrectRejection(1:lastinx)==1) / (sum(FalseAlarm(1:lastinx)==1|CorrectRejection(1:lastinx)==1));

%% Reaction time

HitRT = nanmean(GoRT);
FaRT = nanmean(NoGoRT);

Output = [GoPerf NoGoPerf HitRT FaRT]
% openvar('Output')