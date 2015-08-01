function airpuff_notthirsty
%AIRPUFF_NOTTHIRSTY   Hit rate.
%   AIRPUFF_NOTTHIRSTY compares hit rate in a late session versus hit rate
%   when not water deprived. 50 dB tones are used.
%
%   Se also AIRPUFF_NONOISE, AIRPUNISH2 and AIRPUNISH2B.

% Sessions
sessions = {'n074' '150517b'; ...
    'n075' '150517b'; ...
    'n076' '150517b'; ...
    'n079' '150516b'};

% sessions2 = {'n074' '141126a'; ...
%     'n075' '141126a'; ...
%     'n076' '141127a'; ...
%     'n079' '141206a'};

sessions2 = {'n074' '150114a'; ...   % last
    'n075' '150102a'; ...   % last
    'n076' '150111a'; ...   % last
    'n079' '150114a'};   % last-2

[ehitrate1 efarate1 ehits1 efas1 egos1 enogos1] = main(sessions);   % first session
[ehitrate2 efarate2 ehits2 efas2 egos2 enogos2] = main(sessions2);   % first session

% Chi-square test
for iS = 1:length(sessions)
    [h p] = b_chi2test2([ehits1(iS) egos1(iS)-ehits1(iS)],[ehits2(iS) egos2(iS)-ehits2(iS)],0.05)   % chi square test
end

% Plot
figure;
MN = median([ehitrate1; ehitrate2],2);
SE = [se_of_median(ehitrate1) se_of_median(ehitrate2)];
bar(1:2,MN,'BarWidth',0.7,'FaceColor','none',...
    'EdgeColor',[0.8 0 0])
hold on
E = errorbar(1:2,MN,SE,'r+');
errorbar_tick(E,0)
setmyplot_balazs

keyboard

% -------------------------------------------------------------------------
function [ehitrate efarate ehits efas egos enogos] = main(sessions)

% False alarm rate
NumSessions = size(sessions,1);
[ehitrate, efarate, ehits, efas, egos, enogos] = deal(nan(1,NumSessions));
for iS = 1:NumSessions
    try
        fnm = fullfile('F:\NB_cellbase\',sessions{iS,1},sessions{iS,2},'TE.mat');
        TE = load(fnm);
    catch
        fnm = fullfile('G:\HDB_cellbase\',sessions{iS,1},sessions{iS,2},'TE.mat');
        TE = load(fnm);
    end
    
    NumTrials = length(TE.TrialStart);
    ehits(iS) = sum(TE.Hit==1&TE.StimulusDuration==50);   % number of hits in the engaged period
    efas(iS) = sum(TE.FalseAlarm==1&TE.StimulusDuration==50);   % number of false alarms in the engaged period
    egos(iS) = sum((TE.Hit==1&TE.StimulusDuration==50)|TE.Miss==1&TE.StimulusDuration==50);   % number of go tones in the engaged period
    enogos(iS) = sum((TE.FalseAlarm==1&TE.StimulusDuration==50)|TE.CorrectRejection==1&TE.StimulusDuration==50);   % number of no-go tones in the engaged period
    ehitrate(iS) = ehits(iS) / egos(iS);   % hit rate
    efarate(iS) = efas(iS) / enogos(iS);   % false alarm rate
    
%     figure
%     plot(hitrate,'g')
%     hold on
%     plot(farate,'r')
end

% -------------------------------------------------------------------------
function [h,p,stat] = b_chi2test2(n_d1,n_d2,alpha)
% for binned data

if nargin < 3
    alpha = 0.05;
end

ld1 = sum(n_d1);
ld2 = sum(n_d2);

stat = ld1 * ld2 * sum((n_d1/ld1-n_d2/ld2).^2./(n_d1+n_d2));
q = chi2cdf(stat,length(n_d1)-1);
if q > (1 - alpha)
    h = 1;
else
    h = 0;
end
p = 1 - q;