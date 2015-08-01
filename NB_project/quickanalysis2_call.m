%% Choose cellbase

animalID = 76;
san = num2str(animalID);

if ismember(animalID,[70 78 79])
    cb = 'HDB';
elseif ismember(animalID,[71 72 73 74 75 76 77])
    cb = 'NB';
end


%% sessionspec = behavior, recording and stimulation

choosecb(cb)
quickanalysis2(animalID,'150517b',[1 0 0],'feedbackdelay')

%% waterpuff

choosecb(cb)
quickanalysis2(animalID,'150510a',[1 0 0],'waterpuff')

%% pavlovian

choosecb(cb)
quickanalysis_pavlovian(animalID,'150113a',[1 1 1])

%% APECS 

% sessionspec = behavior, recording and stimulation
% cbroot = getpref('cellbase','datapath');
% sessions = listtag('allsessions');
% NumSessions = size(sessions,1);
% for iS = 1:NumSessions
%     quickanalysis2_APECS(sessions{iS,1},sessions{iS,2})
% end

% quickanalysis2_APECS('cs001','140714a')

%% optic tract

% quickanalysis_optictract2(23,'120105b',0)
% quickanalysis_optictract2(23,'120104b',0)
% quickanalysis_optictract2(23,'120104c',0)
% quickanalysis_optictract2(23,'120102d',0)

%% human

% psych2trialevents_susattn('c:\Balazs\_data\human_susattn\cellbase\ld\130705a\ld_susattn_130705a.mat')
% auditory_gonogo_psychplot3('ld','130705a')