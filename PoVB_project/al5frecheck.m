%% Subrutine for finding optimal starting point for segments of layer 5
%% cortical neurons (constrained in LFP frequency)

freq = [];
next = 1;
for k = 1:200:len-seglen
    ind1 = k;
    ind2 = ind1 + seglen -1;
    vd = vdisc(vdisc>ind1&vdisc<ind2) - ind1;      % localize
    loceeg = eeg(ind1:ind2);
    lfeeg = feeg((ind1-1)/const+1:ind2/const);
    lahee = ahee((ind1-1)/const+1:ind2/const);

    eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
    vdisc2 = round(vd/const);

    cyclen = eegfre(lfeeg,lahee,dsr);    % filter EEG, Hilbert-transform
    freq(next) = 1 / cyclen * 1000;
    next = next + 1;
end

figure
plot(freq)