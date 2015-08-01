animalID = 'n015';
animalID2 = 'nb015';
sessionID = '110623a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

TE1 = solo2trialevents2_auditory_gonogo([fullpth 'data_@auditory_gonogo_balazs_' animalID2 '_' sessionID '.mat']);

TE = load([fullpth 'TrialEvents.mat']);

TE.ITIDistribution = TE1.ITIDistribution;
TE.StimulusDuration = TE1.StimulusDuration;
TE.SoundIntensity = TE1.SoundIntensity;

save([fullpth 'TrialEvents.mat'],'-struct','TE');