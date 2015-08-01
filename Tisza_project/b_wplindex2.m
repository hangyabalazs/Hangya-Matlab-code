function b_wplindex2(wave1, wave2)
%WPLINDEX2   Phase locking indeces.
%   WPLINDEX2(W1,W2) calculates phase coherence, entropy index and mutual 
%   information index for W1 and W2 wavelets.
%
%   See also WCROSSWAVELET.

% Crosswavelet
wave_cross = wave1 .* conj(wave2);

% Phase locking indeces
N = 365;    % 1 year window
T = size(wave_cross,2);T = 2000;
F = size(wave_cross,1);
pow = abs(wave_cross) .^ 2;
phas = angle(wave_cross);
phas1 = angle(wave1);
phas2 = angle(wave2);
bno = fix(exp(0.626+0.4*log(N-1)));   % number of bins
Hmax = log(bno);
PhaseCoh = zeros(F,T-N+1);  % phase coherence
EntrInd = zeros(F,T-N+1);   % entropy index
MutInfInd = zeros(F,T-N+1);   % mutual information index
for t = N:T
    PhaseCoh(1:F,t-N+1) = abs(sum(exp(1).^(i.*phas(1:F,t-N+1:t)),2)./N) .^ 2;  % phase coherence
    
    p = hist(phas(1:F,t-N+1:t)',bno) ./ N;
    warning off
    lp = log(p');
    lp(find(lp<-1000000)) = 0;
    warning backtrace
    EntrInd(1:F,t-N+1) = -sum(p'.*lp,2);
    
    p1 = hist(phas1(1:F,t-N+1:t)',bno) ./ N;
    p2 = hist(phas2(1:F,t-N+1:t)',bno) ./ N;
    p12 = zeros(F,bno,bno);
    for f=1:F
        p1t = (phas1(f,t-N+1:t) + pi) / ((2 * pi) / bno);
        p1tt = fix(p1t) + 1;
        p2t = (phas2(f,t-N+1:t) + pi) / ((2 * pi) / bno);
        p2tt = fix(p2t) + 1;
        for s = 1:length(p1tt)
            p12(f,p1tt(s),p2tt(s)) = (p12(f,p1tt(s),p2tt(s)) + 1) ./ N;
        end
    end
    warning off
    lp1 = log(p1');
    lp1(find(lp1<-1000000)) = 0;
    lp2 = log(p2');
    lp2(find(lp2<-1000000)) = 0;
    lp12 = log(p12);
    lp12(find(lp12<-1000000)) = 0;
    warning backtrace
    p1_entropy = -sum(p1'.*lp1,2);
    p2_entropy = -sum(p2'.*lp2,2);
    joint_entropy = -sum(sum(p12.*lp12,3),2);
    mutual_information(1:F,t-N+1) = p1_entropy + p2_entropy - joint_entropy;
end
EntrInd = (Hmax - EntrInd) ./ Hmax;   % entropy index
MutInfInd = (Hmax - mutual_information) ./ Hmax;   % mutual information index

% Plot
figure
imagesc([366 T],[1 F],PhaseCoh)
figure
imagesc([366 T],[1 F],EntrInd)
figure
imagesc([366 T],[1 F],MutInfInd)