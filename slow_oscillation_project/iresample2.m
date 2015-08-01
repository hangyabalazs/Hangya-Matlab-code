function iresample2(fname)
%IRESAMPLE   Resample multichannel EEG on 1000 Hz.
%   IRESAMPLE(FNAME) resamples multichannel EEG data on 1000 Hz. Input data
%   in FNAME should be sampled on 256 Hz.
%
%   See also DATATRANSFORM and RESAMPLE.

% Load
load(fname)
chno = size(data,2);

% Resample
p = 125;
q = 32;
data2 = zeros(ceil(size(data,1)*p/q),chno);
for k = 1:chno
    chold = data(:,k);
    chnew = resample(chold,p,q);
    data2(:,k) = chnew;
end
data = data2(1:60005,:);

% Save
[pth fn ext] = fileparts(fname);
fname2 = fullfile(pth,[fn '_rs',ext]);
save(fname2,'data')



