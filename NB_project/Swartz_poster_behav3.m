animalID = 'n013';
animalID2 = 'nb013';
sessionID = '110629a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);

[PerfGo1 PerfNoGo1 TimeGo1 TimeNoGo1] = auditory_gonogo_psychplot2(animalID,sessionID,[],impulsivity_filter(TE));
Dprime1 = PerfGo1 - PerfNoGo1;

%%

animalID = 'n013';
animalID2 = 'nb013';
sessionID = '110630a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);

lightoff_trial = find(TE.LightStimulation2==0);
[PerfGo2 PerfNoGo2 TimeGo2 TimeNoGo2] = auditory_gonogo_psychplot2(animalID,sessionID,[],intersect(lightoff_trial,impulsivity_filter(TE)));
Dprime2 = PerfGo2 - PerfNoGo2;

%%

animalID = 'n013';
animalID2 = 'nb013';
sessionID = '110701a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);

lightoff_trial = find(TE.LightStimulation2==0);
[PerfGo3 PerfNoGo3 TimeGo3 TimeNoGo3] = auditory_gonogo_psychplot2(animalID,sessionID,[],intersect(lightoff_trial,impulsivity_filter(TE)));
Dprime3 = PerfGo3 - PerfNoGo3;

%%

animalID = 'n013';
animalID2 = 'nb013';
sessionID = '110704a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);

lightoff_trial = find(TE.LightStimulation2==0);
[PerfGo4 PerfNoGo4 TimeGo4 TimeNoGo4] = auditory_gonogo_psychplot2(animalID,sessionID,[],intersect(lightoff_trial,impulsivity_filter(TE)));
Dprime4 = PerfGo4 - PerfNoGo4;

%%

animalID = 'n013';
animalID2 = 'nb013';
sessionID = '110705a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
TE = load([fullpth 'TrialEvents.mat']);

lightoff_trial = find(TE.LightStimulation2==0);
[PerfGo5 PerfNoGo5 TimeGo5 TimeNoGo5] = auditory_gonogo_psychplot2(animalID,sessionID,[],lightoff_trial);
Dprime5 = PerfGo5 - PerfNoGo5;

%% mean for nb013

PerfGo_nb013 = mean([PerfGo1(1:4); PerfGo2(1:4); PerfGo3(1:4); PerfGo5(1:4)]);
PerfGoSE_nb013 = std([PerfGo1(1:4); PerfGo2(1:4); PerfGo3(1:4); PerfGo5(1:4)]) / 2;

PerfNoGo_nb013 = mean([PerfNoGo1(1:4); PerfNoGo2(1:4); PerfNoGo3(1:4); PerfNoGo5(1:4)]);
PerfNoGoSE_nb013 = std([PerfNoGo1(1:4); PerfNoGo2(1:4); PerfNoGo3(1:4); PerfNoGo5(1:4)]) / 2;

Dprime_nb013 = mean([Dprime1(1:4); Dprime2(1:4); Dprime3(1:4); Dprime5(1:4)]);
DprimeSE_nb013 = std([Dprime1(1:4); Dprime2(1:4); Dprime3(1:4); Dprime5(1:4)]) / 2;

figure
% plot([20 30 40 50],PerfGo_nb013,'r')
errorbar([20 30 40 50],PerfGo_nb013,PerfGoSE_nb013,'g')
hold on
errorbar([20 30 40 50],PerfNoGo_nb013,PerfNoGoSE_nb013,'r')

figure
errorbar([20 30 40 50],Dprime_nb013,DprimeSE_nb013,'g')

TimeGo_nb013 = mean([TimeGo1(1:4); TimeGo2(1:4); TimeGo3(1:4); TimeGo5(1:4)]);
TimeGoSE_nb013 = std([TimeGo1(1:4); TimeGo2(1:4); TimeGo3(1:4); TimeGo5(1:4)]) / 2;

TimeNoGo_nb013 = mean([TimeNoGo1(1:4); TimeNoGo2(1:4); TimeNoGo3(1:4); TimeNoGo5(1:4)]);
TimeNoGoSE_nb013 = std([TimeNoGo1(1:4); TimeNoGo2(1:4); TimeNoGo3(1:4); TimeNoGo5(1:4)]) / 2;

figure
% plot([20 30 40 50],TimeGo_nb013,'r')
errorbar([20 30 40 50],TimeGo_nb013,TimeGoSE_nb013,'g')
hold on
errorbar([20 30 40 50],TimeNoGo_nb013,TimeNoGoSE_nb013,'r')



%%

animalID = 'n014';
animalID2 = 'nb014';
sessionID = '110628a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

[PerfGo1 PerfNoGo1] = auditory_gonogo_psychplot2(animalID,sessionID);

%%

animalID = 'n014';
animalID2 = 'nb014';
sessionID = '110630a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

[PerfGo1 PerfNoGo1] = auditory_gonogo_psychplot2(animalID,sessionID);

%%

animalID = 'n014';
animalID2 = 'nb014';
sessionID = '110628a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

[PerfGo1 PerfNoGo1] = auditory_gonogo_psychplot2(animalID,sessionID);

%%

animalID = 'n015';
animalID2 = 'nb015';
sessionID = '110704a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

[PerfGo1 PerfNoGo1] = auditory_gonogo_psychplot2(animalID,sessionID);

%%

animalID = 'n015';
animalID2 = 'nb015';
sessionID = '110705a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

[PerfGo2 PerfNoGo2] = auditory_gonogo_psychplot2(animalID,sessionID);