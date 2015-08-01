function bgresample(fname)
%BGRESAMPLE   Resample 1250 Hz EEG on 1000 Hz.
%   BGRESAMPLE(FNAME) resamples EEG data on 1000 Hz. Input data in FNAME
%   should be sampled on 1250 Hz.
%
%   See also DATATRANSFORM and RESAMPLE.

% Load
load(fname)
data = Eeg;
chno = size(data,1);

% Resample
p = 4;
q = 5;
data2 = zeros(chno,ceil(size(data,2)*p/q));
for k = 1:chno
    chold = data(k,:);
    chnew = resample(chold,p,q);
    data2(k,:) = chnew;
end
data = data2;

% Save
[pth fn ext] = fileparts(fname);
fname2 = fullfile(pth,[fn '_rs',ext]);
save(fname2,'data')