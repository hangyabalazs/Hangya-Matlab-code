function contrastrun(delta)
%CONTRASTRUN    Runs CONTRASTING on a sequence of files.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code. The figure and the input argument of 
%   CONTRASTING_FOR_CONTRASTRUN (delta) are saved in the result directory during the 
%   execution.
%
%   CONTRASTRUN(DELTA) uses DELTA as input argument of CONTRASTING_FOR_CONTRASTRUN.
%
%   CONTRASTRUN, by itself, uses DELTA = 200 default value.
%
%   See also CONTRASTING_FOR_CONTRASTRUN and CONTRASTING.

% Input arguments check
error(nargchk(0,1,nargin));
if nargin == 0,
    delta = 200;
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
cd([DATAPATH,'relfreq_figs_proba\']);  %Here are the results
npt = input...
    ('Discard existing content or append data while writing contrasting_1.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('contrasting_1.txt','w');
elseif npt == 'a',
    fid = fopen('contrasting_1.txt','a');
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
    
% Computing the input variables
    close all;
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
            
% Contrasting
            [hh] = b_contrasting_for_contrastrun(vdisc,unit,delta);     %contrasting
            
% Saving the plot
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
fprintf(fid,'%s',num2str(delta));
fclose(fid);
close all;
cd(mm);