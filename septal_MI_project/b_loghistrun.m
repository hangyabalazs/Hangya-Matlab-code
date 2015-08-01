function b_loghistrun(m,k,nmax)
%LOGHISTRUN    Runs LOGHIST on a sequence of files.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code.
%
%   LOGHIST(M,K,NMAX) computes logaritmic interspike interval histogram where M and K are
%   constants used for computing the bin borders: M determines the base of the logaritm
%   and K determines the center of the first bin; NMAX is the number of the bins.
%
%   You can specify the input parameters in the program code as well, using LOGHISTRUN by 
%   itself.
%
%   With LOGHISTRUN(A,B) you can give any of the three input arguments setting A to 1, 2
%   or 3 for M, K and NMAX. 
%   Note, that using this form A must be a scalar equal to 1, 2 or 3!
%
%   LOGHISTRUN(M) is LOGHISTRUN with M-based logaritm.
%
%   The histogram is saved in the result directory during the execution. The input parameters
%   M, K, and NMAX are saved in a text file.
%
%   See also LOGHIST and CONTRASTLOGHIST.

% Input arguments check
error(nargchk(0,3,nargin));
switch nargin
case 0
    m = 30;
    k = 0.25;
    nmax = 50;
case 1
    k = 0.25;
    nmax = 50;
case 2
    if m == 1,
        m = k;
        k = 0.25;
        nmax = 50;
    elseif m == 2,
        m = 30;
        nmax = 50;
    elseif m == 3,
        nmax = k;
        m = 30;
        k = 0.25;
    else error('The first input argument must be scalar equal to 1, 2 or 3!');
    end;
end;


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
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
cd([DATAPATH,'loghist_figs_proba\']);  %Here are the results
npt = input...
    ('Discard existing content or append data while writing loghist_proba.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('loghist_proba.txt','w');
elseif npt == 'a',
    fid = fopen('loghist_proba.txt','a');
else
    error('Unexpected answer for input.')
end;

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
        meret=size(data,1);
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
    for p = 1:3,
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
           
% Loghist
            [q,h] = b_loghist_for_loghistrun(m,k,isi,nmax);
            
% Creating the histogram
            hh = figure;
            subplot(2,1,1);
            plot(time,unit,'m');
            subplot(2,1,2);
            bar(q,h);
            
% Saving the histogram
            pont = findstr(filename,'.');
            filenam = filename(1:pont(1)-1);
            switch p
            case 1
                title(['NO THETA ',filenam(1:3),' ',filenam(5:6)]);
                eval(['saveas(hh,''NO_THETA_',filenam,'.fig'')']);
            case 2
                title(['SP THETA ',filenam(1:3),' ',filenam(5:6)]);
                eval(['saveas(hh,''SP_THETA_',filenam,'.fig'')']);
            case 3
                title(['EV THETA ',filenam(1:3),' ',filenam(5:6)]);
                eval(['saveas(hh,''EV_THETA_',filenam,'.fig'')']);
            end;
        end;
    end;
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator

% Saving the input arguments
fprf1 = ['m = ',num2str(m)];
fprf2 = ['k = ',num2str(k)];
fprf3 = ['nmax = ',num2str(nmax)];

fprintf(fid,'%s\n%s\n%s',fprf1,fprf2,fprf3);
fclose(fid);
close all;
cd(mm)