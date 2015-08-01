function imicontrol_snr_calibrate
%IMICONTROL_SNR_CALIBRATE   Arteficial control for shifted mutual information calculation.
%   An interesting try...
%
%   See also IMICONTROL1 and IMISHIFT.

lw = 800;
gs = gausswin(lw) * 10;           % 0.8 s gaussian 'wave'
data = zeros(1,10005);           % 10 second arteficial data, sampled on 1000 Hz
% load('X:\In_Vivo\balazs\_analysis\Ulbert\control_SNR_okosan\lowpropdata.mat')
% data = (dd(1:10005)' - min(dd(1:10005))) / (max(dd(1:10005)) - min(dd(1:10005))) * 10;
data(1,2001:2000+lw) = gs;        % wave at 2 sec. on channel1
data(1,4551:4550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(1,6551:6550+lw) = gs; 
data = data';
gns = [0:0.1:1 2:100 110:10:500];      % adding Gaussian noise
SNR = zeros(1,length(gns));
for k = 1:length(gns)
    pSNR = zeros(1,50);
    for t = 1:50
        data1 = data + rand(size(data)) * gns(k) / 10;
        [FFT,w] = b_fft2(data1,1000,8192);
        pSNR(t) = max(FFT(w>=0.5&w<=4)) / mean(FFT(w>=0.5&w<=40));   % max delta power to all power (Hurtado et al 2004)
    end
    SNR(k) = mean(pSNR);
end
figure
plot(gns,SNR)


1
% save('X:\In_Vivo\raw_data\human_SO\oiti37_lukacs\grid\mat_EEG_12\arteficial_control_snr.mat','data') % save