function b_maxfrun
%MAXFRUN    Runs MAXF on a sequence of files.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   The function creates two figures: a histogram of the top 50 values of the instant
%   frequency and the instant frequency plot. It saves three computed parameters in a text
%   file: the mean, the standard deviation and the median of the top 50 instant frequency 
%   values.
%
%   See also MAXF.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'analysenow5\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH,'maxf_figs\']);  %Here are the results
npt = input...
    ('Discard existing content or append data while writing maxf.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('maxf.txt','w');
elseif npt == 'a',
    fid = fopen('maxf.txt','a');
else
    error('Unexpected answer for input.')
end;

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
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
    
% Computing the input variables
    dt = 0.0001;
    mintafr = 1 / dt;
    unit = data(:,2);
    unit = unit';
    eeg = data(:,1);
    eeg = eeg';
    time = [0:length(unit)-1] * dt; 
    xlimit = [min(time),max(time)];
    figure
    plot(time,unit)
    datinx1 = input('start: ');
    datinx2 = input('end: ');
%     datinx1 = 1;
%     datinx2 = length(unit);
    global IN
    IN = cell(1,12);
    IN{1} = data;
    IN{2} = eeg;
    IN{3} = fname;
    IN{4} = where1;
    IN{5} = datinx1;
    IN{6} = datinx2;
    IN{7} = time;
    IN{8} = unit;
    IN{9} = dt;
    IN{10} = meret;
    IN{11} = mintafr;
    IN{12} = xlimit;
    
% Discrimination
    b_disc;
    global DISC
    id = DISC{1};
    output = DISC{2};
    vdisc = DISC{3};
    kuszob = DISC{4};
    instfrek = DISC{5};
    isi = DISC{6};
    
% MAXF
    [highfr,instfr] = b_maxf_for_maxfrun;
    bins = [50:10:500];
    [n m]=hist(highfr,bins);
    ave_hf = mean(highfr);
    std_hf = std(highfr);
    med_hf = median(highfr);
    max_hf = max(highfr);
    
% Plotting
    h1 = figure;
    bar(m,n)
    h2 = figure;
    plot(instfr)
    
% Saving the plots
    pont = findstr(fname,'.');
    filenam = fname(1:pont(1)-1);
    set(0,'CurrentFigure',h1);
    title([filenam(1:3),' ',filenam(5:6),'   Inst. Fr. Top50 Hist.']);
    eval(['saveas(h1,''MAXF_hist_',filenam(1:6),'.fig'')']);
    set(0,'CurrentFigure',h2);
    title([filenam(1:3),' ',filenam(5:6),'   Instant Frequency']);
    eval(['saveas(h2,''MAXF_instfr_',filenam(1:6),'.fig'')']);
    fprintf(fid,'%s %s %s %s %s\n',filenam(1:6),ave_hf,std_hf,med_hf,max_hf);
    
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
fclose(fid);
close all;
cd(mm);