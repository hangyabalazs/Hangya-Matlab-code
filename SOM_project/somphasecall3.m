function somphasecall3
%SOMPHASECALL3   Phase calculations.
%   SOMPHASECALL3 is functionally similar to SOMPHASECALL but uses an EEG
%   resampling approach instead of interpolation, which results in faster
%   execution. It also calculates Spike Triggered Averages and
%   autocorrelations.
%
%   See also SOMPHASECALL2.

% Directories
global DATAPATH
inpdir = '\\science\Kepecs\Hyun\For Balazs\';
resdir = [DATAPATH '\SOM\Phase_temp\'];
mm = pwd;
cd(resdir)

% Phase histogram
dr = dir(inpdir);
for k = 4:length(dr)
    
    if (dr(k).isdir && strncmp(dr(k).name,'.',1)==0)
        
        % Load EEG
        inpadd = dr(k).name;
        drr = [inpdir inpadd];
        [fl0 fl] = filelist(drr);
        fst = find(strncmp(fl,'CSC',3));
        eegname = fl{fst(1)};
        eegnm = eegname(1:end-4);
        load([inpdir inpadd filesep eegname])
        load([inpdir inpadd filesep 'EVENTS'])  %load events
        eval(['frs = ' eegnm '_SampleFrequencies;'])    % sampling rate
        eval(['eegts = ' eegnm '_TimeStamps;'])    % EEG time stamps
        eval(['clear ' eegnm '_TimeStamps'])    % EEG time stamps
        if max(frs) > 1020
            continue
        end
        
        % Remove LaserStim part
        eval(['eeg0 = ' eegnm '_Samples;']);  % EEG
        eval(['clear ' eegnm '_Samples']);  % EEG
        sr = 1010;
        dt = 1 / sr;
        eeg0 = eeg0(:);
        time = (0:length(eeg0)-1) * dt + eegts(1)/1000000;     % EEG time stamps
        t0 = eegts(1)/1000000;
        if Events_TimeStamps(1) > 1e8
            Events_TimeStamps = Events_TimeStamps / 1000000;
        end
        idx = find(strcmp(Events_EventStrings, 'Cheetah 32 DIO board 0 input TTL (0x0200)'),1);
        idx_start = find(time>Events_TimeStamps(idx),1,'first');
        eeg = eeg0(1:idx_start-1);
        clear eeg0 Events_TimeStamps
        
        % Resample
        eeg = resample(eeg,99,100);
        sr = 1000;
        dt = 1 / sr;
%         time = (0:length(eeg)-1) * dt + eegts(1) / 1000000;     % EEG time stamps
        
        
        %         %plot spectrogram
        %         spectrogram1(eeg,dt);
        %         fname = ['Spectrogram_' inpadd '_' eegnm '.fig'];
        %         saveas(gcf,fname)
        
        % Load unit
        fst0 = find(strncmp(fl,'TT',2));
        for cells = 2:2
            unitname = fl{fst0(cells)};
            unitnm = unitname(1:end-4);
            load([inpdir inpadd filesep unitname])
            TS = TS - eegts(1)/1000000;
            TS = TS(TS<(length(eeg)/sr));
            vdisc = TS' * sr;
            
            % Phase calculation
            [hang hmvl ftm bahee ahee Hhist Hctrl] = somphase(eeg,vdisc,sr,3,10);   % theta
            fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(Hhist,fnm)
%             fnm = ['THETA_PHASEctrl_' inpadd '_' eegnm '_' unitnm '.fig'];
%             saveas(Hctrl,fnm)
            if ~isempty(bahee)
                [Z,p_rayleigh,U,p_rao] = b_rao(bahee');
            end
            disp(p_rao)
            %fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            %save(fnm,'hang','hmvl','ftm','bahee')
            
            betaL = 20;     % beta band boundaries
            betaU = 35;
            [hang, hmvl, ftm, bahee ahee Hhist Hctrl] = somphase(eeg,vdisc,sr,betaL,betaU);   % beta
            %cd('beta')
            fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(Hhist,fnm)
            %fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            %save(fnm,'hang','hmvl','ftm','bahee')
            %cd ..
            clear ahee
            
            % Spike Triggered Average
            wn = 2 * sr;    % 2 sec. window
            [StaTheta StaIndexTheta1 StaIndexTheta2 nn] = astanorm(round(vdisc),eeg,wn);
            titlestr = [inpadd ' ' eegnm ' ' unitnm];
            titlestr(titlestr=='_') = ' ';
            H = stafig(StaTheta,StaIndexTheta1,StaIndexTheta2,nn,wn,sr,titlestr);
            fnm = ['THETA_STA_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(H,fnm);
            
            % Autocorrelation
            wn = sr;
            if ~isempty(TS)
                [H1 H2 trsc] = somccg(TS,TS,wn);  % autocorr.
                figure(H1)
                title(titlestr)
                fnm = ['THETA_AUTO_' inpadd '_' eegnm '_' unitnm '.fig'];
                saveas(H1,fnm);
            end
            
            % Wavelet
            dsr = 500;
            const = sr / dsr;
            seglen = 10 * 60 * dsr;        % 10 min. long segments
            len0 = length(eeg);
            lenr = floor(len0/seglen);
            ind1 = 1:seglen:len0;
            ind2 = ind1 + seglen -1;
            for es = 1:lenr
                eegl = eeg(ind1(es):ind2(es));
                vdiscl = TS(TS*sr>ind1(es)&TS*sr<ind2(es));
                eegw = eegl(1:const:end);     % downsampling EEG
                vdiscw = round(vdiscl*dsr);   % downsampling unit
                len = length(eegw);
                [wave_eeg,f] = eegwavelet(eegw,dsr);        % EEG wavelet
                im_eeg = 10 * log10(wave_eeg);      % dB scale
                dt = 1 / dsr;
                wavetime = (0:len-1) * dt;
                clear wave_eeg
                
                Hw = figure;
                S1 = subplot(3,1,1);
                imagesc(im_eeg)
                set(gca,'CLim',[0 20])
                mx = max(im_eeg(:));
                x_lim = get(gca,'XLim');
                y_lim = get(gca,'YLim');
                xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
                ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
                text(xcoord,ycoord,num2str(mx),'Color',[1 1 1]);
                ff = round(f*100) / 100;
                time = round(wavetime*100) / 100;
                b_rescaleaxis('Y',ff)
                b_rescaleaxis('X',time)
                title(titlestr)
                setappdata(S1,'scalex',wavetime)
                setappdata(S1,'scaley',f(f>0.5))
                b_zoomset_for_wavelet
                S2 = subplot(3,1,2);
                dneeg = (eegw - mean(eegw)) / std(eegw);
                etime = 1:size(im_eeg,2);
                plot(etime,dneeg)
                xlim([etime(1) etime(end)])
                hold on
                flt = fir1(1024,[betaL betaU]/(dsr/2),'band');      % filtering in beta range
                fdneeg = filtfilt(flt,1,dneeg);
                plot(etime,fdneeg,'r')
                S3 = subplot(3,1,3);
                dnunit = zeros(1,etime(end));
                dnunit(vdiscw) = 1;
                plot(etime,dnunit)
                xlim([etime(1) etime(end)])
                setappdata(S1,'subplothandle',[S2 S3])
                fns = ['WAVELET_' inpadd '_' eegnm '_' unitnm '_' num2str(es) '.fig'];  % save
                saveas(Hw,fns)
            end
            
            close all
        end
    end
end
cd(mm)



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [hang hmvl ftm bahee ahee H Hc] = somphase(eeg,vdisc,sr,fr1,fr2)
%SOMPHASE   Phase histogram.
%   [A R FTM ANGS AHEE] = SOMPHASE(EEG,VD,SRFR1,FR2) calculates phase
%   histograms for neuronal spikes (VD) relative to EEG sampled at SR.
%   Phase is calculated via Hilbert-transform of the EEG filtered between
%   FR1 and FR2 Hz. Output arguments are mean phase (A), mean vector length
%   (R), first trigonometric moment (FTM), all unit phase values (ANGS) and
%   all EEG phase values (AHEE).
%
%   [A R FTM ANGS AHEE H HC] = SOMPHASE(EEG,VD,SRFR1,FR2) returns figure
%   handles of unit phase histogram and a control phase histogram, from all
%   EEG phase values.
%
%   See also CZPHASE2.

% Filtering EEG
nqf = sr / 2;
flt = fir1(2048,[fr1 fr2]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);
figure
plot((eeg-mean(eeg))/std(eeg))
hold on
plot(feeg,'g')

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 1/fr2 s
% 3. discard cicles longer then 1/fr1 s
% fn = find(-diff(ahee)>2*pi-0.3);
fn0 = valuecrossing(1:length(ahee),ahee',0,'down');
fn = round(fn0);
sd = std(feeg);
inx = [];
vdisc(vdisc<0|vdisc>length(eeg)) = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx find(vdisc>fn0(k)&vdisc<fn0(k+1))];
        plot(fn(k)+1:fn(k+1),seeg,'r')
    end
end
vdisc(inx) = [];

% Phase angles - Hilbert
bahee = ahee(round(vdisc));
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Plot phase histogram
angs = bahee;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm,xout] = histc(angs*180/pi,edges);   % phase histogram
nm = nm(1:end-1);
H = figure;
B = bar(cnts,nm'/length(angs));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(length(angs))];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

[nm,xout] = histc(ahee*180/pi,edges);   % control phase histogram (all EEG phase values)
nm = nm(1:end-1);
Hc = figure;
B = bar(cnts,nm'/length(ahee));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])

% -------------------------------------------------------------------------
function H = stafig(sta,sta_index1,sta_index2,nn,wn,sr,titlestr)

time = linspace(-wn/sr/2,wn/sr/2,length(sta));
H = figure;
plot(time,sta,'LineWidth',1.5)
ach = allchild(H);     % figure title
ax = findobj(ach,'type','axes');
title(ax(end),titlestr)
x_lim = xlim;
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

% -------------------------------------------------------------------------
function [wave,f] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 0;
disp('No padding!')
dj = 0.02;    
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
[wave,period,scale,coi] = b_wavelet_pow3(dat,dt,pad,dj,s0,j1,mother,param,mis);

% -------------------------------------------------------------------------
function b_zoomset_for_wavelet
%ZOOMSET_FOR_WAVELET    Assigns ZOOM_FOR_WAVELET as figure ButtonDownFcn.
%
%   See also ZOOM2 and ZOOMSET.

% Get handles
ax = gca;
plt = findobj(ax,'Type','line');
ptch = findobj(ax,'Type','patch');
im = findobj(ax,'Type','image');

% Set application data
if ~isempty(plt)
    x_data = get(plt,'XData');
    y_data = get(plt,'YData');
    setappdata(ax,'x_data',x_data)
    setappdata(ax,'y_data',y_data)
end
x_lim = get(ax,'XLim');
y_lim = get(ax,'YLim');
setappdata(ax,'x_lim',x_lim)
setappdata(ax,'y_lim',y_lim)

% Set ButtonDownFcn
set(ax,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
if ~isempty(plt)
    set(plt,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
elseif ~isempty(ptch)
    set(ptch,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
elseif ~isempty(im)
    set(im,'ButtonDownFcn','b_zoom_for_wavelet_linkaxes')
end