function b_ivsirun
%IVSIRUN    Runs IVSI on a sequence of files.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code. The function creates two kind of fig files: the IVSI
%   figure and another figure showing the bursts found by BURST_CLS.
%
%   See also IVSI, BURST_CLS and WAVEPHASE.

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
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
cd([DATAPATH,'ica_figs_proba\']);  %Here are the results
wts = input('Do you wish to save the number of clusters? /yes: Enter, no: 0/','s');

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    wbo =o;     %o for waitbar
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
            
% Ivsi            
            [burstoc,intervals,minimums,o,H] = b_ivsi_main_for_ivsirun(wts);
            
% Creating the plots            
            pont = findstr(filename,'.');
            filenam = filename(1:pont(1)-1);
            switch p
            case 1
                if isempty(o) == 0,
                    set(0,'CurrentFigure',o);
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),'   IVSI']);
                    eval(['saveas(o,''NO_THETA_IVSI_',filenam,'.fig'')']);
                end;
                if isempty(H) == 0,
                    set(0,'CurrentFigure',H);
                    title(['NO THETA ',filenam(1:3),' ',filenam(5:6)]);
                    eval(['saveas(H,''NO_THETA_bursts_',filenam,'.fig'')']);
                end;
            case 2
                if isempty(o) == 0,
                    set(0,'CurrentFigure',o);
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6),'   IVSI']);
                    eval(['saveas(o,''SP_THETA_IVSI_',filenam,'.fig'')']);
                end;
                if isempty(H) == 0,
                    set(0,'CurrentFigure',H);
                    title(['SP THETA ',filenam(1:3),' ',filenam(5:6)]);
                    eval(['saveas(H,''SP_THETA_bursts_',filenam,'.fig'')']);
                end;
            case 3
                if isempty(o) == 0,
                    set(0,'CurrentFigure',o);
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6),'   IVSI']);
                    eval(['saveas(o,''EV_THETA_IVSI_',filenam,'.fig'')']);
                end;
                if isempty(H) == 0,
                    set(0,'CurrentFigure',H);
                    title(['EV THETA ',filenam(1:3),' ',filenam(5:6)]);
                    eval(['saveas(H,''EV_THETA_bursts_',filenam,'.fig'')']);
                end;
            end;
        end;
    end;
    waitbar(wbo/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
close all;
cd(mm)