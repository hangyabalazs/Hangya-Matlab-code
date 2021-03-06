function imicontrol_snr4
%IMICONTROL_SNR4   Arteficial control for shifted mutual information calculation.
%   IMICONTROL_SNR4 generates arteficial control data by shifting a gaussian
%   window. A variable amount of Gaussian noise is added to the channels to
%   mimic the effect of differences in signal-tonoise ratio. See the
%   comments in the code for the precise structure of the
%   data generated!
%
%   See also IMICONTROL1 and IMISHIFT.

% Load number of propagating waves andnoise levels
load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\gauss4.mat')

rtim = round(randTIM/6);
csnr = randNoise;

% Generate data
lw = 800;
gs = gausswin(lw) * 29.3395 * 4.4410;           % 0.8 s gaussian 'wave' (amplitude adjusted to real amplitude ratios)
sd = 10004;
data = zeros(20,sd);           % 10 second arteficial data, sampled on 1000 Hz
data(1,8101:8100+lw) = gs;        % wave at 2 sec. on channel1
data(1,101:100+lw) = gs;        % wave at 2 sec. on channel1
data(2,5701:5700+lw) = gs;        % wave at 2 sec. on channel1
data(2,851:850+lw) = gs;        % wave at 2 sec. on channel1
data(3,2801:2800+lw) = gs;        % ch4 -> ch2 & ch5 (800 ms)
data(4,3101:3100+lw) = gs;
data(5,3501:3500+lw) = gs;        % ch2 & ch6 -> ch7 (900 ms)
data(6,1851:1850+lw) = gs;       % ch7 -> ch12 (850 ms)
data(7,1551:1550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(8,4851:4850+lw) = gs;        % ch7 -> ch8 (850 ms)
data(8,6551:6550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(10,4251:4250+lw) = gs;        % ch7 -> ch8 (850 ms)
data(11,3951:3950+lw) = gs;        % ch7 -> ch8 (850 ms)
data(12,351:350+lw) = gs;        % ch7 -> ch8 (850 ms)
data(12,4551:4550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(12,5401:5400+lw) = gs;        % ch7 -> ch8 (850 ms)
data(13,8701:8700+lw) = gs;       % ch10 -> ch20 (600 ms)
data(14,1001:1000+lw) = gs;       % ch19 -> ch15 (700 ms)
data(15,2101:2100+lw) = gs;       % ch20 -> ch10 (600 ms)
data(15,6801:6800+lw) = gs;       % ch19 -> ch15 (700 ms)
data(17,7101:7100+lw) = gs;       % ch12 -> ch20 (750 ms)
data(17,5101:5100+lw) = gs;       % ch12 -> ch20 (750 ms)
data(17,9101:9100+lw) = gs;       % ch12 -> ch20 (750 ms)
data(18,5951:5950+lw) = gs;       % ch12 -> ch20 (750 ms)
data(18,7501:7500+lw) = gs;       % ch15 -> ch20 (700 ms)
data(19,2401:2400+lw) = gs;       % ch20 -> ch19 (800 ms)
data(19,6201:6200+lw) = gs;       % ch20 -> ch19 (800 ms)
data(20,1301:1300+lw) = gs;       % ch12 -> ch20 (750 ms)

% Add white noise levels
data = data';
% lline = length(find(w>4));
% lline2 = length(find(w<=4));
lline2 = 41;
lline = 10004 - lline2;
for k = 1:20
    ns = csnr(k);
    ns2 = ns / lline;
    ns3 = [zeros(1,lline2) ones(1,lline)*ns2];
    ng = rand(1,length(ns3)) * 2 * pi - pi;
    cf = sqrt(ns3) .* exp(i*ng);
    cf = cf';
    leny = length(cf);
    if isequal(mod(leny,2),0)
        hlf = leny / 2;
        cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
    else
        hlf = (leny + 1) / 2;
        cf(hlf+1:end)=flipud(conj(cf(2:hlf)));
    end
    ceeg = ifft(cf);
    ceeg = real(ceeg);
    data(:,k) = data(:,k) + ceeg;
end
1
save('X:\In_Vivo\raw_data\human_SO\oiti37_lukacs\grid\mat_EEG_12\arteficial_control_snr3.mat','data') % save