function bggammamod_rob(data,Laps)
%BGGAMMAMOD_ROB   Modulation of gamma power by theta phase.
%   Analysis based on the article of Robert Komorowski (Tort AB, Komorowski
%   RW, Manns JR, Kopell NJ, Eichenbaum H. (2009) Theta-gamma coupling
%   increases during the learning of item-context associations. Proc Natl
%   Acad Sci U S A. 2009 Nov 23.) See BGGAMMAMOD_RAYLEIGH2 for details.
%
%   See also BGGAMMAMOD_RAYLEIGH2.

% Results directory
global DATAPATH
resdir = [DATAPATH 'Wheel2\Gammamod_rayleighbox\'];
mm = pwd;
cd(resdir)

% Get trial EEG
lap_start = find(Laps.StartLaps)';
lap_end = [lap_start(2:end)-1 size(Laps.TrialType,1)];
lap_mazesection = round(resample(Laps.MazeSection,4,5));
lap_center = round((lap_start+lap_end)/2);
lap_type = Laps.TrialType(lap_center);
lap_behav = Laps.BehavType(lap_center);
lap_start = round(lap_start*4/5);  % indices are transformed to match new samp. rate
lap_end = [lap_start(2:end)-1 size(data,2)];
lap_center = round((lap_start+lap_end)/2);
lt = length(lap_start);
% lt = 1;
% error_WE = [];
% futureL_WE = [];
% futureR_WE = [];
dsglR = 1:20:18000;
dsglL = 1:20:18000;
dsghR = 1:20:18000;
dsghL = 1:20:18000;
aMIglR = [];
aMIglL = [];
aMIghR = [];
aMIghL = [];
aMIglR_err = [];
aMIglL_err = [];
aMIghR_err = [];
aMIghL_err = [];
aMIglR_fL = [];
aMIglL_fL = [];
aMIghR_fL = [];
aMIghL_fL = [];
aMIglR_fR = [];
aMIglL_fR = [];
aMIghR_fR = [];
aMIghL_fR = [];
modix = zeros(lt,10) * NaN;
modix2 = zeros(lt,10) * NaN;
modix3 = zeros(lt,10) * NaN;
modix4 = zeros(lt,10) * NaN;
tpixR = zeros(lt,10) * NaN;
lgpixR = zeros(lt,10) * NaN;
hgpixR = zeros(lt,10) * NaN;
tpixC = zeros(lt,10) * NaN;
lgpixC = zeros(lt,10) * NaN;
hgpixC = zeros(lt,10) * NaN;
for k = 1:lt
    k
    if ~isequal(lap_behav(k),1) || isequal(lap_type(k),0)
        continue
    end
    tS = lap_start(k);
    tE = lap_end(k);
%     time = linspace(tS,tE,length(trials{k}.Distance)*4/5+1);
%     d1 = trials{k}.Distance(1);
%     d2 = trials{k}.Distance(end-1);
%     pd = linterp(linspace(tS,tE,length(trials{k}.Distance)),...
%         trials{k}.Distance,time(2:end-1));
%     distance = [d1 pd d2];   % distance from the start point
    eegR = data(10,tS:tE);   % shank: 4 or 10
    eegL = data(4,tS:tE);
    maze_section = lap_mazesection(tS:tE);
    
% Filter
    sr = 1000;
    nqf = sr / 2;
    flt = fir1(1024,[4 12]/nqf,'band');      % bandpass filtering: theta frequency band bounderies: 8-10 Hz
    feeg_thetaR = filtfilt(flt,1,eegR);
    feeg_thetaR = (feeg_thetaR - mean(feeg_thetaR)) / std(feeg_thetaR);
    feeg_thetaL = filtfilt(flt,1,eegL);
    feeg_thetaL = (feeg_thetaL - mean(feeg_thetaL)) / std(feeg_thetaL);
    hilb_thphaseR = angle(hilbert(feeg_thetaR));
    hilb_thphaseL = angle(hilbert(feeg_thetaL));
    flt = fir1(1024,[30 40]/nqf,'band');      % bandpass filtering: low gamma frequency band bounderies: 30-40 Hz
    feeg_lowR = filtfilt(flt,1,eegR);
    feeg_lowR = (feeg_lowR - mean(feeg_lowR)) / std(feeg_lowR);
    feeg_lowL = filtfilt(flt,1,eegL);
    feeg_lowL = (feeg_lowL - mean(feeg_lowL)) / std(feeg_lowL);
    hilb_glabsR = abs(hilbert(feeg_lowR));
    hilb_glabsL = abs(hilbert(feeg_lowL));
    flt = fir1(1024,[80 120]/nqf,'band');      % bandpass filtering: high gamma frequency band bounderies: 80-120 Hz
    feeg_highR = filtfilt(flt,1,eegR);
    feeg_highR = (feeg_highR - mean(feeg_highR)) / std(feeg_highR);
    feeg_highL = filtfilt(flt,1,eegL);
    feeg_highL = (feeg_highL - mean(feeg_highL)) / std(feeg_highL);
    hilb_ghabsR = abs(hilbert(feeg_highR));
    hilb_ghabsL = abs(hilbert(feeg_highL));
    
% Wavelet
%     [waveR f] = eeg_wavelet(eegR,sr);   % EEG sampled on 1000 Hz
%     [waveL f] = eeg_wavelet(eegL,sr);
%     wave_cross = waveR .* conj(waveL);
%     waveaR = abs(waveR) .^ 2;
%     waveaL = abs(waveL) .^ 2;
%     wavea_cross = abs(wave_cross) .^ 2;
%     clear waveR waveL wave_cross
%     pwind1t = find(f>4,1,'last');    % frequency band bounderies
%     pwind2t = find(f<12,1,'first');
%     pwind1lg = find(f>30,1,'last');
%     pwind2lg = find(f<40,1,'first');
%     pwind1hg = find(f>80,1,'last');
%     pwind2hg = find(f<120,1,'first');
%     tpowerR = mean(waveaR(pwind2t:pwind1t,:));
%     tpowerL = mean(waveaL(pwind2t:pwind1t,:));
%     tpowerC = mean(wavea_cross(pwind2t:pwind1t,:));
%     lgpowerR = mean(waveaR(pwind2lg:pwind1lg,:));
%     lgpowerL = mean(waveaL(pwind2lg:pwind1lg,:));
%     lgpowerC = mean(wavea_cross(pwind2lg:pwind1lg,:));
%     hgpowerR = mean(waveaR(pwind2hg:pwind1hg,:));
%     hgpowerL = mean(waveaL(pwind2hg:pwind1hg,:));
%     hgpowerC = mean(wavea_cross(pwind2hg:pwind1hg,:));
        
% Modulation vs. distance
    thetaR = nrm1(hilb_thphaseR); % normalization
    thetaL = nrm1(hilb_thphaseL);
    gammalR = nrm1(hilb_glabsR);
    gammalL = nrm1(hilb_glabsL);
    gammahR = nrm1(hilb_ghabsR);
    gammahL = nrm1(hilb_ghabsL);

% Rayleigh-test
    [MItRglR MItRglL MItLglR MItLglL MItRghR MItRghL MItLghR MItLghL pglR pglL pghR pghL] = ...
        lrayleigh(thetaR,thetaL,gammalR,gammalL,gammahR,gammahL,1000,10);
% Plot
    tm = linspace(0,tE-tS,length(MItRglR));
    tmm = linspace(0,tE-tS,length(maze_section));
    for mseg = 1:10
%         inx = find(maze_section==mseg);
%         tpixR(k,mseg) = nanmedian(tpowerR(inx));
%         lgpixR(k,mseg) = nanmedian(lgpowerR(inx));
%         hgpixR(k,mseg) = nanmedian(hgpowerR(inx));
%         tpixC(k,mseg) = nanmedian(tpowerC(inx));
%         lgpixC(k,mseg) = nanmedian(lgpowerC(inx));
%         hgpixC(k,mseg) = nanmedian(hgpowerC(inx));
        
        
        cfdf = MItRglR;
%         cfdf2 = cfdf;
%         cfdf2(isnan(pglR)|pglR>0.01) = NaN;
        inx = find(maze_section==mseg);
        if isempty(inx)
            continue
        end
        sne = seglim(inx);
        inx2 = [];
        for inc = 1:size(sne,2)
            pinx = find(tm>tmm(sne(1,inc))&tm<tmm(sne(2,inc)));
            inx2 = [inx2 pinx];
        end
        mcf = cfdf(inx2);
%         mcf2 = mcf(mcf>0);
%         mcf3 = cfdf2(inx2);
%         mcf4 = mcf3(mcf3>0);
        modix(k,mseg) = nanmedian(mcf);
%         modix2(k,mseg) = nanmedian(mcf2);
%         modix3(k,mseg) = nanmedian(mcf3);
%         modix4(k,mseg) = nanmedian(mcf4);
    end
    
    
   
% Save
%     switch lap_type(k)
%         case {3,4}
%             cd error
% 
%         case 2
%             cd futureL
% 
%         case 1
%             cd futureR
% 
%         otherwise
%             keyboard
%     end
%     fn = ['g01_m13_t' num2str(k) '_MIglRt.fig'];
%     saveas(HglRt,fn)
%     fn = ['g01_m13_t' num2str(k) '_MIglLt.fig'];
%     saveas(HglLt,fn)
%     fn = ['g01_m13_t' num2str(k) '_MIghRt.fig'];
%     saveas(HghRt,fn)
%     fn = ['g01_m13_t' num2str(k) '_MIghLt.fig'];
%     saveas(HghLt,fn)
    
    cd ..
end
% save vars aMIglR aMIglL aMIghR aMIghL aMIglR_err aMIglL_err aMIghR_err aMIghL_err...
%     aMIglR_fL aMIglL_fL aMIghR_fL aMIghL_fL aMIglR_fR aMIglL_fR aMIghR_fR aMIghL_fR
cd(mm)



% -------------------------------------------------------------------------
% RAYLEIGH'S TEST
% -------------------------------------------------------------------------
function [MItRglR MItRglL MItLglR MItLglL MItRghR MItRghL MItLghR MItLghL ...
    pglR pglL pghR pghL] = ...
    lrayleigh(thetaR,thetaL,gammalR,gammalL,gammahR,gammahL,WindowSize,Overlap)

[k1 k2] = size(thetaR);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);

% Gamma thresholds
TglR = mean(gammalR) + 2*std(gammalR);
TglL = mean(gammalL) + 2*std(gammalL);
TghR = mean(gammahR) + 2*std(gammahR);
TghL = mean(gammahL) + 2*std(gammahL);

% Rayleigh's test
ovlp = Overlap;
MItRglR = zeros(1,maxi*ovlp-ovlp+1);
MItRglL = zeros(1,maxi*ovlp-ovlp+1);
MItLglR = zeros(1,maxi*ovlp-ovlp+1);
MItLglL = zeros(1,maxi*ovlp-ovlp+1);
MItRghR = zeros(1,maxi*ovlp-ovlp+1);
MItRghL = zeros(1,maxi*ovlp-ovlp+1);
MItLghR = zeros(1,maxi*ovlp-ovlp+1);
MItLghL = zeros(1,maxi*ovlp-ovlp+1);
pglR = zeros(1,maxi*ovlp-ovlp+1) * NaN;
pglL = zeros(1,maxi*ovlp-ovlp+1) * NaN;
pghR = zeros(1,maxi*ovlp-ovlp+1) * NaN;
pghL = zeros(1,maxi*ovlp-ovlp+1) * NaN;
for i = 1:maxi*ovlp-ovlp+1        % ABS LOOP
    inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlapping windows!
    inx1 = round(inx1);
    inx2 = inx1 + winlen - 1;
    
    yj1 = thetaR(inx1:inx2);
    yj2 = gammalR(inx1:inx2);
    pbang = find(yj2>TglR);
    if length(pbang) < 10
        MItRglR(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tRglR = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItRglR(i) = z;
        else
           MItRglR(i) = NaN;
        end
    end
    
    yj1 = thetaR(inx1:inx2);
    yj2 = gammalL(inx1:inx2);
    pbang = find(yj2>TglL);
    if length(pbang) < 10
        MItRglL(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tRglL = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItRglL(i) = z;
        else
            MItRglL(i) = NaN;
        end
    end
    
    yj1 = thetaL(inx1:inx2);
    yj2 = gammalR(inx1:inx2);
    pbang = find(yj2>TglR);
    if length(pbang) < 10
        MItLglR(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tLglR = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItLglR(i) = z;
        else
            MItLglR(i) = NaN;
        end
        [u2 pw] = b_watsontwo(bang_tRglR',bang_tLglR');
        pglR(i) = pw(2);
    end
    
    yj1 = thetaL(inx1:inx2);
    yj2 = gammalL(inx1:inx2);
    pbang = find(yj2>TglL);
    if length(pbang) < 10
        MItLglL(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tLglL = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItLglL(i) = z;
        else
            MItLglL(i) = NaN;
        end
        [u2 pw] = b_watsontwo(bang_tRglL',bang_tLglL');
        pglL(i) = pw(2);
    end
    
    yj1 = thetaR(inx1:inx2);
    yj2 = gammahR(inx1:inx2);
    pbang = find(yj2>TghR);
    if length(pbang) < 10
        MItRghR(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tRghR = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItRghR(i) = z;
        else
            MItRghR(i) = NaN;
        end
    end
    
    yj1 = thetaR(inx1:inx2);
    yj2 = gammahL(inx1:inx2);
    pbang = find(yj2>TghL);
    if length(pbang) < 10
        MItRghL(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tRghL = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItRghL(i) = z;
        else
            MItRghL(i) = NaN;
        end
    end
    
    yj1 = thetaL(inx1:inx2);
    yj2 = gammahR(inx1:inx2);
    pbang = find(yj2>TghR);
    if length(pbang) < 10
        MItLghR(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tLghR = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItLghR(i) = z;
        else
            MItLghR(i) = NaN;
        end
        [u2 pw] = b_watsontwo(bang_tRghR',bang_tLghR');
        pghR(i) = pw(2);
    end
    
    yj1 = thetaL(inx1:inx2);
    yj2 = gammahL(inx1:inx2);
    pbang = find(yj2>TghL);
    if length(pbang) < 10
        MItLghL(i) = NaN;
    else
        bang = yj1(pbang);
        bang_tLghL = bang;
        [z p] = rayt(bang);
        if p < 0.001
            MItLghL(i) = z;
        else
            MItLghL(i) = NaN;
        end
        [u2 pw] = b_watsontwo(bang_tRghL',bang_tLghL');
        pghL(i) = pw(2);
    end
    
end     % end of abs loop

% -------------------------------------------------------------------------
function [z p] = rayt(bang)

n = length(bang);
ftm = sum(exp(1).^(j*bang)) / n;    % first trigonometric moment
mrl = abs(ftm);     % mean resultant length
z = n * (mrl ^ 2);  % Rayleigh's Z statistic
p = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
    (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));



% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function [wave,f] = eeg_wavelet(dat,sr)
%EEGWAVELET   Wavelet calculation.

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 1;
dj = 0.08;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);



% -------------------------------------------------------------------------
% NORMALIZATION
% -------------------------------------------------------------------------
function N = nrm1(pN)

N = pN - min(pN);
N = N / max(N);



% -------------------------------------------------------------------------
% INTERPOLATION
% -------------------------------------------------------------------------
function yi = llinterp(x,y,xi)
%LINTERP   1D linear interpolation.
%   YI = LINTERP(X,Y,XI) interpolates Y = f(XI), where f is given in (X,Y).
%
%   See also INTERP1Q.

% Interpolation
yi = zeros(1,length(xi));
for k = 1:length(xi)    % linear interpolation
    inx1 = find(x<xi(k),1,'last');
    if isempty(inx1) || inx1==length(x)
        continue
    end
    inx2 = inx1 + 1;
    yi(k) = y(inx1) + (y(inx2) - y(inx1)) * (xi(k)-x(inx1)) / (x(inx2) - x(inx1));
end



% -------------------------------------------------------------------------
% SEGMENTATION
% -------------------------------------------------------------------------
function sne = seglim(ip)

% Makes disjuct intervals from a set of points
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
prepa = zeros(2,lenfdr+1);
prepa(1,1) = ip(1);
for t = 1:lenfdr
    prepa(2,t) = ip(fdr(t));
    prepa(1,t+1) = ip(fdr(t)+1);
end
prepa(2,end) = ip(end);
sne = prepa;