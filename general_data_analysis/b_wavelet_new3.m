%WAVELET  1D Wavelet transform with optional singificance testing
%
%   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM,MIS)
%
%   Computes the wavelet transform of the vector Y (length N),
%   with sampling rate DT.
%
%   By default, the Morlet wavelet (k0=6) is used.
%   The wavelet basis is normalized to have total energy=1 at all scales.
%
%
% INPUTS:
%
%    Y = the time series of length N.
%    DT = amount of time between each Y value, i.e. the sampling time.
%
% OUTPUTS:
%
%    WAVE is the WAVELET transform of Y. This is a complex array
%    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%    The WAVELET power spectrum is ABS(WAVE)^2.
%    Its units are sigma^2 (the time series variance).
%
%
% OPTIONAL INPUTS:
% 
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects (see COI below).
%
%    DJ = the spacing between discrete scales. Default is 0.25.
%         A smaller # will give better scale resolution, but be slower to plot.
%
%    S0 = the smallest scale of the wavelet.  Default is 2*DT.
%
%    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%    MOTHER = the mother wavelet function.
%             The choices are 'MORLET', 'PAUL', or 'DOG'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%            For 'PAUL' this is m (order), default is 4.
%            For 'DOG' this is m (m-th derivative), default is 2.
%
%    MIS = the maximal intresting scale.
%        Wavelet is computed only at the scales of interest in order to reduce
%        memory usage and execution time. Default value is J1+1.
%
%
% OPTIONAL OUTPUTS:
%
%    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
%           to the SCALEs.
%
%    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
%            where J1+1 is the total # of scales.
%
%    COI = if specified, then return the Cone-of-Influence, which is a vector
%        of N points that contains the maximum period of useful information
%        at that particular time.
%        Periods greater than this are subject to edge effects.
%        This can be used to plot COI lines on a contour plot by doing:
%
%              contour(time,log(period),log(power))
%              plot(time,log(coi),'k')
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%----------------------------------------------------------------------------
function [wave,period,scale,coi] = b_wavelet_new3(Y,dt,pad,dj,s0,J1,mother,param,mis);
%WAVELET_NEW3   Wavelet optimalized for memory usage.
%   Some memory conserving modifications were made on the original WAVELET function.
%   To reach the best memory perfofmance, use WAVELET_NEW2 (execution will be slower).
%
%   See also WAVELET_NEW2, WAVELET_POW2 and WAVELET_POW3.

% Input argument check
if (nargin < 9), mis = -1; end
if (nargin < 8), param = -1;, end
if (nargin < 7), mother = -1;, end
if (nargin < 6), J1 = -1;, end
if (nargin < 5), s0 = -1;, end
if (nargin < 4), dj = -1;, end
if (nargin < 3), pad = 0;, end
if (nargin < 2)
	error('Must input a vector Y and sampling time DT')
end

n1 = length(Y);

if (s0 == -1), s0 = 2 * dt;, end
if (dj == -1), dj = 1 ./ 4.;, end
if (J1 == -1), J1 = fix((log(n1*dt/s0)/log(2))/dj);, end
if (mother == -1), mother = 'MORLET';, end

% Construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2)+0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

% Construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];
k = k .* ((2 .* pi) / (n * dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

% Compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

% Construct SCALE array & empty PERIOD & WAVE arrays
scale = s0 * 2 .^ ((0:J1) * dj);
period = scale;

% Free memory
clear Y pad dj s0 x base2 n expr
clear functions

% Loop through all scales and compute transform
if mis == -1
    mis = J1 + 1;
end
wave = zeros(mis,n1);  % define the wavelet array
for a1 = 1:mis
	[daughter,fourier_factor,coi,dofmin] = wave_bases(mother,k,scale(a1),param);	
	ifd = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
    ifd = ifd(1:n1);
    wave(a1,:) = ifd;
end

period = fourier_factor*scale;
coi = coi * dt * [1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]