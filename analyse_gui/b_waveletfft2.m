function b_waveletfft2(wave,f,samprate)
%WAVELETFFT   Calculates Fourier transform of every scale of wavelet.
%   Three input arguments are needed for this function: the wavelet transformated
%   eeg or unit, the scale - frequency matrix and the sampling rate (the new
%   sampling rate, if resampled).
%   The result is shown in three images for three frequency bands (under theta,
%   theta, over theta), treating the values as hights above a plane.
%
%   See also FFT, FFT2, WAVELET, WAVELET_NEW2 and WAVELET_NEW3.

%Input arguments check
error(nargchk(3,3,nargin));

% Find theta band
fnd = find(f>6);
pwind1 = fnd(end);
fnd = find(f<3);
pwind2 = fnd(1);

% Power and axis
pow = abs(wave).^2;
[sw1 sw2] = size(pow);
X = 'X';
Y = 'Y';

% Calculate and plot fft
figure

S1 = subplot(3,1,1);  % 1st subplot
w1 = pow(1:pwind1,:);
fr1 = f(1:pwind1);
p1 = zeros(pwind1,sw2);
f1 = zeros(1,sw2);
for i = 1:pwind1
    [p1(i,:) f1(1,:)] = b_fft2(w1(i,:),samprate,sw2,2);
end
fnd = find(f1<15);
lim1 = fnd(end);
pet1 = p1(:,4:lim1);
axes(S1)
imagesc(pet1);
b_rescaleaxis(X,f1)     % rescaling
b_rescaleaxis(Y,fr1)
setappdata(S1,'scalex',f1);     % assign zoom as buttondown function
setappdata(S1,'scaley',fr1);
b_zoomset_for_wavelet

S2 = subplot(3,1,2);    % 2nd subplot
w2 = pow(pwind1+1:pwind2,:);
fr2 = f(pwind1+1:pwind2);
sz = pwind2-pwind1;
p2 = zeros(sz,sw2);
f2 = zeros(1,sw2);
for i = 1:sz
    [p2(i,:) f2(1,:)] = b_fft2(w2(i,:),samprate,sw2,2);
end
fnd = find(f1<2);
lim2 = fnd(end);
pet2 = p2(:,4:lim2);
axes(S2)
imagesc(pet2);
b_rescaleaxis(X,f2)     % rescaling
b_rescaleaxis(Y,fr2)
setappdata(S2,'scalex',f2);     % assign zoom as buttondown function
setappdata(S2,'scaley',fr2);
b_zoomset_for_wavelet

S3 = subplot(3,1,3);    % 3rd subplot
w3 = pow(pwind2+1:end,:);
fr3 = f(pwind2+1:end);
sz = sw1-pwind2;
p3 = zeros(sz,sw2);
f3 = zeros(1,sw2);
for i = 1:sz
    [p3(i,:) f3(1,:)] = b_fft2(w3(i,:),samprate,sw2,2);
end
fnd = find(f1<1);
lim3 = fnd(end);
pet3 = p3(:,4:lim3);
axes(S3);
imagesc(pet3);
b_rescaleaxis(X,f3)     % rescaling
b_rescaleaxis(Y,fr3)
setappdata(S3,'scalex',f3);     % assign zoom as buttondown function
setappdata(S3,'scaley',fr3);
b_zoomset_for_wavelet