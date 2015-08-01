function xf = bandpassFilter(x,Fs,Fp1,Fp2)
% function xf = bandpassFilter(X,FS,FP1,FP2)
%
% Bandpass filter for the signal x.  An acausal fft 
% algorithm is applied (i.e. no phase shift). The filter functions is         
% constructed from a Hamming window (default window used in "fir2" Matlab function). 
% to avoid ripples in the frequency reponse (windowing is a smoothing in frequency domain)
%
% Fs : sampling frequency
%
% The passbands (Fp1 Fp2) frequencies are defined in Hz as
%                  ----------                      
%                /|         | \
%               / |         |  \
%              /  |         |   \
%             /   |         |    \
%   ----------    |         |     ----------------- 
%                 |         |
%           Fs1  Fp1       Fp2   Fs2            
%
% DEFAULTS values
% Fs1 = Fp1 - 0.5 in Hz
% Fs2 = Fp2 + 0.5 in Hz
%
%
% If NO OUTPUTS arguments are assigned the filter function H(f) and
% impulse response are plotted. 
%
% NOTE: for long data traces the filter is very slow.
%
% EXEMPLE 
%    x= sin(2*pi*12*[0:1/200:10])+sin(2*pi*30*[0:1/200:10])
%    y=bandpassFilter(x,200,5,20);    bandpass filter between 5 and 20 Hz
%------------------------------------------------------------------------
% Originally produced by the Helsinki University of Technology,
% Adapted by Mariecito SCHMUCKEN 2001
%------------------------------------------------------------------------

%Default values in Hz
Fs1 = Fp1 - 0.5; 
Fs2 = Fp2 + 0.5;

if size(x,1) == 1
    x = x';
end
% Make x EVEN
Norig = size(x,1); 
if rem(Norig,2)
    x = [x' zeros(size(x,2),1)]';                
end

% Normalize frequencies  
Ns1 = Fs1/(Fs/2);
Ns2 = Fs2/(Fs/2);
Np1 = Fp1/(Fs/2);
Np2 = Fp2/(Fs/2);

% Construct the filter function H(f)
N = size(x,1);
Nh = N/2;

B = fir2(N-1,[0 Ns1 Np1 Np2 Ns2 1],[0 0 1 1 0 0]); 
% Make zero-phase filter function
H = abs(fft(B));  
IPR = real(ifft(H));

if size(x,2) > 1
    for k=1:size(x,2)
        xf(:,k) = real(ifft(fft(x(:,k)) .* H'));
    end
    xf = xf(1:Norig,:);
else
    xf = real(ifft(fft(x') .* H));
    xf = xf(1:Norig);
end
x=x(1:Norig); 

% if NO OUTPUT argument then plots
if nargout == 0 
    f = Fs*(0:Nh-1)/(N);
    freqz(IPR,1,f,Fs);    
    figure, subplot(2,1,1)
    plot(f,H(1:Nh));
    xlim([0 2*Fs2])
    title('Filter function H(f)')
    xlabel('Frequency (Hz)')
    subplot(2,1,2)
    plot((1:Nh)/Fs,IPR(1:Nh))
    xlim([0 2/Fp1])
    xlabel('Time (sec)')
    ylim([min(IPR) max(IPR)])
    title('Impulse response')
    figure, subplot(211),
    periodogram(x,hamming(Norig),1024,Fs);
    subplot(212), periodogram(xf,hamming(Norig),1024,Fs);
end

return