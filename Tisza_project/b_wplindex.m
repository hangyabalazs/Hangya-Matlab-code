function b_wplindex(wave_cross)
%WPLINDEX   Phase locking indeces.
%   WPLINDEX(CW) calculates phase coherence and entropy index for CW
%   crosswavelet.
%
%   See also WCROSSWAVELET.

% Phase locking indeces
N = 365;    % 1 year window
T = size(wave_cross,2);T = 2000;
F = size(wave_cross,1);
pow = abs(wave_cross) .^ 2;
phas = angle(wave_cross);
bno = fix(exp(0.626+0.4*log(N-1)));   % number of bins
Hmax = log(bno);
PhaseCoh = zeros(F,T-N+1);  % phase coherence
EntrInd = zeros(F,T-N+1);   % entropy index
for t = N:T
    PhaseCoh(1:F,t-N+1) = abs(sum(exp(1).^(i.*phas(1:F,t-N+1:t)),2)./N) .^ 2;
    
    p = hist(phas(1:F,t-N+1:t)',bno) ./ N;
    warning off
    lp = log(p');
    lp(find(lp<-1000000)) = 0;
    warning backtrace
    EntrInd(1:F,t-N+1) = -sum(p'.*lp,2);
end


EntrInd = (Hmax - EntrInd) ./ Hmax;

% Plot
figure
imagesc([366 T],[1 F],PhaseCoh)
figure
imagesc([366 T],[1 F],EntrInd)