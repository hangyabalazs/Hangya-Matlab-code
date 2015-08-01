function hzinstfreq
%HZINSTFREQ    Instantenous frequency.
%   HZINSTFREQ calculates instantenous frequency for discriminated data 
%   files.
%
%   See also RAPHEPREVIEW.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR2
where = DATADIR2;    %Here are the discriminated data files
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'Data\Instfreq2\']);  %Here are the results

% Progress indicator
wb = waitbar(0,'Running HZINSTFREQ...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Instantenous frequency
sr = 10000;     % sampling rate
for o = 248:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ffnm = [where fname];
    load(ffnm);     % load data
    
    lenu = round(length(eeg)/10);
    dsr = 1000;
    instfrek = zeros(1,lenu);      % instantenous frequency
    vdisc = round(vdisc/10);          % downsample on 1000 Hz
    vdisc = vdisc(vdisc>0);
    isi = diff(vdisc);
    for k = 1:length(vdisc)-1
        instfrek(vdisc(k):vdisc(k+1)) = 1 / isi(k);
    end
    instfrek(1:vdisc(1)-1) = 1 / vdisc(1);
    instfrek(vdisc(end):lenu) = 1 / (lenu - vdisc(end));
    instfrek2 = instfrek * dsr;
    ffr = length(vdisc(vdisc<30*dsr)) / 30;
    sif = smooth(instfrek2,'linear',5001);
    figure
    plot(instfrek2)
    hold on
    plot(sif,'g')
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.2,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,num2str(ffr))
    title(titlestr)
    fn = [fname(1:end-4) '_INSTFREQ'];
    saveas(gcf,fn);
    close all
    waitbar(o/sf)
end
close(wb)