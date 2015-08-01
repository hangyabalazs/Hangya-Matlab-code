function b_spnodist
%SPNODIST   Spike number distribution.

% Input arguments check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
global DATADIR
% where = [DATAPATH,'DATA\analysenow3\'];    %Here are the data files
where = DATADIR;
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'Entropy\Spikecount\']);  %Here are the results
if ~b_isdir2('thres')
    mkdir thres
end

% Progress indicator
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Import
spno = [];
for o = 1:sf    %"CELL CYCLE"
    fname = files(o).name;
    ffnm = [where fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    data = b_load_data(ffnm);
    
    datinx1 = 1;             %first point of the interval
    datinx2 = 1500000;       %last point of the interval
    try
        b_imp(fname,where,data,datinx1,datinx2);
        fname
    catch
        waitbar(o/sf)   %Progress indicator
        disp(['Short registration: ' fname]);
        continue
    end
    
% Thresholding & Discrimination
    [T,segmlen] = b_thres4;
        
    b_disc2(T,segmlen);
    global DISC
    id = DISC{1};
    output = DISC{2};
    vdisc = DISC{3};
        
    cd thres
    str = ['THRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'];
    eval(['save ' str ' T segmlen vdisc']);
    cd ..

% Spike number count
    spno(end+1) = length(vdisc);
    waitbar(o/sf)   %Progress indicator
end
close(wb)

save SpikeNumber spno
hb = fix(exp(0.626+0.4*log(length(spno))));     % number of bins
[x,y] = hist(spno,hb);
B = bar(y,x);
eval(['saveas(B,''SPNODIST.fig'')']);