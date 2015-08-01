%% simulation of Reviewer #2 (expectancy project)


%%%%%%%%%%%%%%%%%%%%%
clear

Fs    = 1000;            % sampling rate in Hz
nyq   = Fs/2;            % nyquist in Hz
band  = [0.5 3.0];       % delta frequency band (Hz)
Wn    = band./nyq;       % normalized frequencies

tvals = [-1499:1500];    % [-1.5s before to 1.5s after]


%%%
% approximate a target-evoked potential
%%%
k = normpdf([-100:100],0,50);  
eeg = zeros(1,3000);
eeg(1550:1600) = +5;
eeg(1750:1800) = -12;
eeg = conv2(eeg,k);
eeg = eeg(49:3048);

[b,a] = fir1(256,Wn,'bandpass');  % bandpass filter

figure(1); clf;
subplot(2,1,1);
hold on;

for i=1:400
  neeg = eeg+40*(rand(1,3000)-0.5); % add noise to each trial
  neeg = filtfilt(b,1,neeg);        % apply bandpass filter
  plot(tvals,neeg);                 % show trial
  h_eeg = hilbert(neeg);            % apply hilbert transform
  p_eeg = angle(h_eeg);             % calculate instant. phase
  t0phase(i) = p_eeg(1500);         % keep track of phase at t=0
end

plot(tvals,eeg,'r');
set(gca,'XLim',[-1000 1000]);
plot([0 0],get(gca,'YLim'),'k');
title('Histogram of phase angles');

subplot(2,2,3);
hist(t0phase,[-pi:0.1:pi]);
set(gca,'XLim',[-pi pi]);
xlabel ('phase (rad)');
title('Histogram of phase angles');
subplot(2,2,4)
rose(t0phase);
title('Rose plot of phase angles');

[Z,p1,U,p2] = b_rao(t0phase)