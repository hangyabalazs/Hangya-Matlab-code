function imicontrol_snr2
%IMICONTROL_SNR2   Arteficial control for shifted mutual information calculation.
%   IMICONTROL_SNR2 generates arteficial control data by shifting a gaussian
%   window. A variable amount of Gaussian noise is added to the channels to
%   mimic the effect of differences in signal-tonoise ratio. See the
%   comments in the code for the precise structure of the
%   data generated!
%
%   See also IMICONTROL1 and IMISHIFT.

load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\rand_TIM_SNR.mat')
load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\calibration_zero_peak.mat')
gns0 = gns;
SNR0 = SNR;
load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\calibration_one_peak.mat')
gns1 = gns;
SNR1 = SNR;
load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\calibration_two_peaks.mat')
gns2 = gns;
SNR2 = SNR;
load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\calibration_three_peaks.mat')
gns3 = gns;
SNR3 = SNR;

rtim = round(randTIM/6);
snrs = zeros(1,20);
for k = 1:20
    eval(['cgns = gns' num2str(rtim(k)) ';']);
    eval(['csnr = SNR' num2str(rtim(k)) ';']);
    inx = find(csnr<randSNR(k),1,'first');
    if isequal(inx,1)
        snrs(k) = 2;
    elseif isempty(inx)
        snrs(k) = 150;
    else
        snrs(k) = gns(inx);
    end
end
snrs(snrs>150) = 150;

lw = 800;
gs = gausswin(lw) * 10;           % 0.8 s gaussian 'wave'
sd = 10004;
data = zeros(20,sd);           % 10 second arteficial data, sampled on 1000 Hz
data(1,2001:2000+lw) = gs;        % wave at 2 sec. on channel1
data(4,2801:2800+lw) = gs;        % ch4 -> ch2 & ch5 (800 ms)
data(5,3101:3100+lw) = gs;
data(6,3701:3700+lw) = gs;        % ch2 & ch6 -> ch7 (900 ms)
data(7,3351:3350+lw) = gs;        % ch7 -> ch8 (850 ms)
data(7,4551:4550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(7,7001:7000+lw) = gs;        % ch7 -> ch8 (850 ms)
data(8,4751:4750+lw) = gs;        % ch7 -> ch8 (850 ms)
data(9,1551:1550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(10,4251:4250+lw) = gs;        % ch7 -> ch8 (850 ms)
data(11,3951:3950+lw) = gs;        % ch7 -> ch8 (850 ms)
data(12,4851:4850+lw) = gs;       % ch7 -> ch12 (850 ms)
data(13,8101:8100+lw) = gs;       % ch20 -> ch10 (600 ms)
data(13,8701:8700+lw) = gs;       % ch10 -> ch20 (600 ms)
data(14,1801:1800+lw) = gs;       % ch19 -> ch15 (700 ms)
data(15,6801:6800+lw) = gs;       % ch19 -> ch15 (700 ms)
data(16,7101:7100+lw) = gs;       % ch12 -> ch20 (750 ms)
data(17,5301:5300+lw) = gs;       % ch12 -> ch20 (750 ms)
data(19,2301:2300+lw) = gs;       % ch20 -> ch19 (800 ms)
data(20,5301:5300+lw) = gs;       % ch12 -> ch20 (750 ms)
data(20,7501:7500+lw) = gs;       % ch15 -> ch20 (700 ms)

data = data';
for k = 1:20
    data(:,k) = data(:,k) + rand(size(data(:,k))) * snrs(k) / 10;
end
1
save('X:\In_Vivo\raw_data\human_SO\oiti37_lukacs\grid\mat_EEG_12\arteficial_control_snr2.mat','data') % save