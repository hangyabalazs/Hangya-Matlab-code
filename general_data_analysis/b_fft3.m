function [f,w] = b_fft3(x,sf,nfft,dim)
%FFT3   Fourier Spectrum
%   F = FFT3(X) calculates complex Fourier spectrum F.
%
%   [F,W] = FFT3(X,SF) calculates Fourier F and frequency vector W using
%   sampling frequency SF.
%
%   [F,W] = FFT3(X,SF,NFFT) calculates NFFT-point Fourier F and frequency
%   vector W.
%
%   [F,W] = FFT3(X,SF,NFFT,DIM) applies the FFT operation across the 
%   dimension DIM.
%
%   The default value of NFFT (if NFFT is empty or not specified) is 
%   size(X,DIM). The default value of DIM (if DIM is empty or not 
%   specified) is 1. If no output arguments are required, FFT3 plots
%   the result (using W frequency vector, if a nonempty SF is given).
%   When X is a row matrix and DIM is not specified, X is converted to
%   column matrix.
%
%   See also FFT and PERIODOGRAM.

% Input arguments check
error(nargchk(1,4,nargin));
if nargout > 1
    if nargin < 2
        warning('If no sampling frequency is given, frequency vector will be empty.')
        w = [];
    end
    if nargin > 1 & isempty(sf)
        warning('If sampling frequency is empty, frequency vector will be empty as well.')
        w = [];
    end
end
if nargin < 4 | (nargin > 3 & isempty(dim))
    dim = 1;
    [s1 s2] = size(x);
    if s1 == 1
        x = x';
    end
end
if nargin < 3 | (nargin > 2 & isempty(nfft))
    nfft = size(x,dim);
end

% Calculate fft
f = fft(x,nfft,dim);

% Generate the frequency vector
if nargin > 1 & ~isempty(sf)
    w = (0 : 1./nfft : 1-1./nfft) .* sf;
end

% If no output arguments required, plot result; clear unrequired output
% arguments
if nargout < 1
    if nargin > 1 & ~isempty(sf)
        plot(w,f)
    else
        plot(f)
    end
    clear f w
end
if nargout == 1
    clear w
end