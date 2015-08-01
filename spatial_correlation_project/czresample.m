%% sampling rate

sr = 1 / CSC6.interval

%% p and q

% p = 132;
% q = 25;
% p = 99;
% q = 50;
p= 12;
q = 25;
sr * p / q


%% resample

eeg = resample(CSC6.values,p,q);

%% save

save('F:\balazs\_analysis\Czurko\czurko_EEG\EEG_CSC6_hux026-day08-tr3-base-01_rs.mat','eeg')

%% clear

clear