function b_clusterrun
%CLUSTERRUN Runs WARD_CLUSTERING on a sequence of files.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code.
%
%   The invoked WARD_CLUSTERING saves a histogram of the intraburst instant frequency and
%   saves also two matrices: burst length in 'bl' and intraburst spike number in 'sn'.
%   Output is stored in d:\matlab6_1\data\burst_clusters\ufos2 directory.
%
%   WARD_CLUSTERING calls a subfunction called PLOTBURST who saves the 'burst plot' in the 
%   result directory.
%
%   See also BURST_CLS and WARD_CLUSTERING.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'analysenow3\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
cd([DATAPATH,'relfreq_figs_proba\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    meret=size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2)==1,
        data=[data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
   
% Searching for the matching datinx file
    cl = fname(1:6);
    filename = [cl,'.txt'];
    pont = findstr(filename,'.');
    filenam = filename(1:pont(1)-1);
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
    ww = 0;
    close all;
    
% Computing the input variables
    p = 3;  %running only on evoked theta interval
    datinx1 = v(2*p-1); %first point of the interval
    datinx2 = v(2*p); %last point of the interval
    
    dt = 0.0001;
    mintafr = 1 / dt;
    if datinx1(1) ~= 0 & datinx2(1) ~= 0,
        ww = ww + 1;
        unit = data(datinx1:datinx2,2);
        unit = unit';
        eeg = data(datinx1:datinx2,1);
        eeg = eeg';
        time = [0:length(unit)-1] * dt; 
        xlimit = [min(time),max(time)];
        global IN
        IN = cell(1,12);
        IN{1} = data;
        IN{2} = eeg;
        IN{3} = fname;
        IN{4} = pathname;
        IN{5} = datinx1;
        IN{6} = datinx2;
        IN{7} = time;
        IN{8} = unit;
        IN{9} = dt;
        IN{10} = meret;
        IN{11} = mintafr;
        IN{12} = xlimit;
        
% Discrimination                
        kuszob = v(6+ww);
        b_disc(kuszob);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};
        kuszob = DISC{4};
        instfrek = DISC{5};
        isi = DISC{6};
        
% Ward's clustering
        b_ward_clustering
        
    end;
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
close all;
cd(mm);