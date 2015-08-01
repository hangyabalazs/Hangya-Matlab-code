function nospect
%NOSPECT   Wavelet spectrum and segment selection.
%   NOSPECT calculates wavelet power spectrum for EEG files of the NO
%   project. Theta segments are selected based on thresholding of
%   point-wise total power. Theta segments are assessed as segments below
%   the mean and of at least half second length. Non-theta segments are
%   established as segments above mean + SD. Transitios are defined as
%   segments between the two thresholds. Non-theta segments with less than
%   half second apart are united.
%
%   See also EEGWAVELET and THETA_MAIN2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
resdir = [DATAPATH 'NO\Wavelet\'];
resdir2 = [resdir 'Segments\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSPECT...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 5000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
newstep = 50;
dsr = sr / newstep;    % downsample on 100 Hz
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname(1:end-4),'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    
    load([inpdir fname])
    eeg = hEEG.values;
    [poweeg,f] = eegwavelet(eeg(1:newstep:end),dsr);    % EEG wavelet
    figure
    imagesc(poweeg)
    set(gca,'CLim',[0 50])
    b_rescaleaxis('Y',f)
    title(titlestr)
    fn = [fname(1:end-4) '_EEGWAVELET.jpg'];
    saveas(gcf,fn);
    
%     pwind1_theta = find(f>6,1,'last');    % theta frequency band bounderies
%     pwind2_theta = find(f<2.5,1,'first');
%     pwind1_delta = find(f>2.5,1,'last');    % delta frequency band bounderies
%     pwind2_delta = find(f<0.5,1,'first');
%     thetapower = mean(poweeg(pwind1_theta:pwind2_theta,:));
%     deltapower = mean(poweeg(pwind1_delta:pwind2_delta,:));
%     tpd = thetapower ./ deltapower;
%     dpt = deltapower ./ thetapower;
    spe = sum(poweeg);      % all power
    spe = spe / std(spe);
    thr = mean(spe) + std(spe);     % define thresholds
    thr2 = mean(spe) + 0.5 * std(spe);
    thr3 = mean(spe);
%     figure
%     sde = std(eeg(1:newstep:end));
%     plot(eeg(1:newstep:end)/sde)
%     hold on
%     plot(spe,'r')
%     eeg2 = eeg(1:5:end)';
%     flt = fir1(4096,[2 8]/500,'band');
%     feeg_theta = filtfilt(flt,1,eeg2)';
%     sdfe = std(feeg_theta(1:10:end));
%     plot(feeg_theta(1:10:end)/sdfe,'c')
%     line([0 length(spe)],[thr thr],'Color','green')
%     line([0 length(spe)],[thr3 thr3],'Color','green','LineStyle',':')
%     line([0 length(spe)],[thr2 thr2],'Color','green','LineStyle',':')
    ip = find(spe<thr3);    % theta segments
    if ~isempty(ip)
        th = segm(ip,newstep,length(eeg));
    end
    th = short_killer(th,sr);   % skip segments shorter than half sec (~2 cycles)
        
    ip = find(spe>thr);     % non-theta segments
    if ~isempty(ip)
        noth1 = segm(ip,newstep,length(eeg));
    end
    
    ip = find(spe<thr&spe>thr3);    % transition segments
    if ~isempty(ip)
        noth2 = segm(ip,newstep,length(eeg));
    end
    
    [noth1, noth2] = changeseg(noth1,noth2,sr);     % unite non-theta segments with shorter than half sec gap
    
    figure      % plot classification result
    plot(eeg)
    hold on
    for k = 1:size(th,2)
        plot(th(1,k):th(2,k),eeg(th(1,k):th(2,k)),'r')
    end
    for k = 1:size(noth1,2)
        plot(noth1(1,k):noth1(2,k),eeg(noth1(1,k):noth1(2,k)),'c')
    end
    for k = 1:size(noth2,2)
        plot(noth2(1,k):noth2(2,k),eeg(noth2(1,k):noth2(2,k)),'g')
    end
    
    fn = [resdir2 fname(1:end-4) '_SEGMENTS.fig'];      % save
    saveas(gcf,fn);
    fn = [resdir2 fname(1:end-4) '_SEGMENTS.mat'];
    ThetaSegments = th;
    NonThetaSegments = noth1;
    Transitions = noth2;
    save(fn,'ThetaSegments','NonThetaSegments','Transitions')
    
    close all
    waitbar(o/sf)
end
close(wb)
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
function [pow,f] = eegwavelet(dat,sr)

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
mif = 0.4;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;

% ----------------------------------------------------------------------------------------
function th = segm(ip,newstep,lenu)
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
preth = zeros(2,lenfdr+1);
preth(1,1) = ip(1);
for t = 1:lenfdr
    preth(2,t) = ip(fdr(t));
    preth(1,t+1) = ip(fdr(t)+1);
end
preth(2,end) = ip(end);
th = preth .* newstep;
th(end) = min(th(end),lenu);

% -------------------------------------------------------------------------
function segments2 = uniteseg(segments,sr)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    if (la-fi)/sr < 0.5
        [fnx fny] = find(segments2==fi);
        segments2(fnx,fny) = segments(2,k+1);
        segments2 = [segments2(1:2,1:fny) segments2(1:2,fny+2:end)];
    end
end

% ----------------------------------------------------------------------------------
function segments = short_killer(segments,sr)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<0.5*sr);         % leaving segments shorter than 0.5 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];

% ----------------------------------------------------------------------------------
function [seg12 seg2] = changeseg(seg1,seg2,sr)

seg12 = seg1;
for k = 1:size(seg1,2)-1
    fi = seg1(2,k);
    la = seg1(1,k+1);
    if la - fi < 0.5 * sr
        pre = find(seg2(1,:)>fi,1,'first');
        if seg2(1,pre) < la
            seg2(:,pre) = [];
        end
        [fnx fny] = find(seg12==fi);
        seg12(fnx,fny) = seg1(2,k+1);
        seg12 = [seg12(1:2,1:fny) seg12(1:2,fny+2:end)];
    end
end