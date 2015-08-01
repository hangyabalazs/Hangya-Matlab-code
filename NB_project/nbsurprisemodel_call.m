function nbsurprisemodel_call
%NBSURPRISEMODEL_CALL   Call NBSURPRISEMODEL_GONOGO_UPDATING.
%   NBSURPRISEMODEL_CALL calls NBSURPRISEMODEL_GONOGO_UPDATING to calculate
%   average performance. See NBSURPRISEMODEL_GONOGO_UPDATING for details.
%
%   See also NBSURPRISEMODEL_GONOGO_UPDATING.

% Run HMM
NumMice = 16;  % number of model mice
NumInt = 4;   % number of model difficulty levels
[GoHitUpdate NoGoHitUpdate GoFAUpdate NoGoFAUpdate] = deal(nan(NumMice,NumInt));
for k = 1:NumMice
    [GoHitUpdate(k,:) NoGoHitUpdate(k,:) GoFAUpdate(k,:) NoGoFAUpdate(k,:)] = ...
        nbsurprisemodel_gonogo_updating;
    close all
end

% Plot average
keyboard
SI = [20 30 40 50];
figure
plot(SI,mean(GoHitUpdate),'g')
hold on
errorbar(SI,mean(GoHitUpdate),nanse(GoHitUpdate),'g')
plot(SI,mean(NoGoHitUpdate),'r')
errorbar(SI,mean(NoGoHitUpdate),nanse(NoGoHitUpdate),'r')

figure
plot(SI,mean(GoFAUpdate),'g')
hold on
errorbar(SI,mean(GoFAUpdate),nanse(GoFAUpdate),'g')
plot(SI,mean(NoGoFAUpdate),'r')
errorbar(SI,mean(NoGoFAUpdate),nanse(NoGoFAUpdate),'r')

keyboard
SI = [0 0.5 1];
GoHitUpdate2 = [GoHitUpdate(:,1) mean([GoHitUpdate(:,2) GoHitUpdate(:,3)],2) GoHitUpdate(:,4)];
NoGoHitUpdate2 = [NoGoHitUpdate(:,1) mean([NoGoHitUpdate(:,2) NoGoHitUpdate(:,3)],2) NoGoHitUpdate(:,4)];
GoFAUpdate2 = [GoFAUpdate(:,1) mean([GoFAUpdate(:,2) GoFAUpdate(:,3)],2) GoFAUpdate(:,4)];
NoGoFAUpdate2 = [NoGoFAUpdate(:,1) mean([NoGoFAUpdate(:,2) NoGoFAUpdate(:,3)],2) NoGoFAUpdate(:,4)];
figure
plot(SI,mean(GoHitUpdate2),'g')
hold on
errorbar(SI,mean(GoHitUpdate2),nanse(GoHitUpdate2),'g')
plot(SI,mean(NoGoHitUpdate2),'r')
errorbar(SI,mean(NoGoHitUpdate2),nanse(NoGoHitUpdate2),'r')

figure
plot(SI,mean(GoFAUpdate2),'g')
hold on
errorbar(SI,mean(GoFAUpdate2),nanse(GoFAUpdate2),'g')
plot(SI,mean(NoGoFAUpdate2),'r')
errorbar(SI,mean(NoGoFAUpdate2),nanse(NoGoFAUpdate2),'r')