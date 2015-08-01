function SNIFF = process_sniffing(Raw_LFP,Fs,Flow)
% processed sniffing trace = process_sniffing(raw_LFP,Sampling frequency,Flow)
% this function will process the sniffing signal
% raw lfp should be denoised and detrended
% It will lowpass filter the data at Flow
% SPR 11-11-2008

if nargin<2,
    Fs = 1000;
end
Nyq = Fs/2;
if nargin<3,
    Flow = 100;
end

[b,a] = butter(5,Flow/Nyq,'low');
SNIFF = filtfilt(b,a,Raw_LFP);