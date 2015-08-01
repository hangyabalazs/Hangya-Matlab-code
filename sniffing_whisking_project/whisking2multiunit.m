function W_TS = whisking2multiunit(W,Fs,freqc,sigma)
%WHISKING2MULTIUNIT   Converts whisker EMG into multi-unit activity.
%
%  W_TS = whisking2multiunit(W,Fs,freqc,sigma)
%
%   Fs - sampling frequency
%   freqc - lowpass cutoff frequency (e.g. 20Hz)
%   sigma - threshold for unit (e.g. 1.5)

% AK 11/2008

%sigma = 1.5;
%freqc = 20;
%Fs=1000;

W = W(:)';
[m,n]=size(W);
time=(1:n)/Fs;

fc=freqc/(Fs/2);
[b,a]=butter(7,fc,'low');

for iT=1:m
    
    Wi_filt=filtfilt(b,a,W(iT,:));
    
    Vi=W(iT,:)-Wi_filt;
    threshold = sigma*std(Vi);
    
    pos=peakdetect(abs(Vi),threshold);
    
    W_TS{iT}=time(pos)';
    
end