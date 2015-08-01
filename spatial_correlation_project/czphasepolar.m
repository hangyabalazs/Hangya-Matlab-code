%% polar plot

load('F:\balazs\_analysis\Czurko\czurko_EEG\phase\resubmission\pos_nocrit\FTMS_pos_rayleigh.mat')
figure
P1 = polar(angle(ftms_rayleigh),abs(ftms_rayleigh),'ro');
set(P1,'MarkerSize',16)
hold on

load('F:\balazs\_analysis\Czurko\czurko_EEG\phase\resubmission\neg_nocrit\FTMS_neg_rayleigh.mat')
P2 = polar(angle(ftms_rayleigh),abs(ftms_rayleigh),'b.');
set(P2,'MarkerSize',20)

%% polar plot for unpaired interneurons

% load('F:\balazs\_analysis\Czurko\czurko_EEG\phase\pos\FTMS_pos_rayleigh.mat')
figure
% P1 = polar(angle(ftms_rayleigh),abs(ftms_rayleigh),'ro');
% set(P1,'MarkerSize',16)
% hold on

load('F:\balazs\_analysis\Czurko\czurko_EEG\phase\resubmission\nopair\FTMS_pos_rayleigh.mat')
tt = 0.3 * exp(pi*i);
P0 = polar(angle(tt),abs(tt),'k.');
hold on
P1 = polar(angle(ftms_rayleigh),abs(ftms_rayleigh),'k.');
set(P1,'MarkerSize',20)
delete(P0)

%% polar plot pyr. cells

load('F:\balazs\_analysis\Czurko\czurko_EEG\phase\resubmission\pyr\FTMS_pos_rayleigh.mat')
figure
P1 = polar(angle(ftms_rayleigh),abs(ftms_rayleigh),'g.');
set(P1,'MarkerSize',16)