function b_relfreqrun
%RELFREQRUN Runs RELFREQ on a sequence of files.
%   This function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code. The function creates two kind of fig files: the RELFREQ
%   figure and Fourier spectrum of it.
%
%   See also FFT and RELFREQ.

%Input arguments check
error(nargchk(0,0,nargin));

%Import
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
% fid = fopen('loghist_3.txt','w');

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = load(ffnm);
    meret = size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2) == 1,
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
   
%Searching for the matching datinx file
    cl = fname(1:6);
    filename = [cl,'.txt'];
    pont = findstr(filename,'.');
    filenam = filename(1:pont(1)-1);
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
    ww = 0;
    close all;
    
%Computing the input variables
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
            
%Discrimination                
            kuszob = v(6+ww);
            b_disc(kuszob);
            global DISC
            id = DISC{1};
            output = DISC{2};
            vdisc = DISC{3};
            kuszob = DISC{4};
            instfrek = DISC{5};
            isi = DISC{6};
            
%Computing relative frequency and fft of relative frequency
            [quociens,dvd] = b_relfreq_for_relfreqrun;
            dvdh = dvd/2;
            loc = vdisc(1:end-1) + dvdh;
            loc = fix(loc);
            val = zeros(1,length(unit));
            val(loc) = quociens;
            pquo = fft(val,65536);
            Pyy = pquo .* conj(pquo) / 65536;
            fquo = 10000 * (0:65536) / 65536;
                        
%Plotting
            close all;
            h1 = figure;
            subplot(2,1,1);
            plot(val);
            subplot(2,1,2);
            plot(unit);
            h2 = figure;
            subplot(2,1,1);
            plot(fquo(1:100),Pyy(1:100));
            subplot(2,1,2);
            plot(unit);
            
%Saving the plots
            pont = findstr(filename,'.');
            filenam = filename(1:pont(1)-1);
            switch p
            case 1
                set(0,'CurrentFigure',h1);
                title(['NO THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
                eval(['saveas(h1,''NO_THETA_quo_',filenam,'.fig'')']);
                set(0,'CurrentFigure',h2);
                title(['NO THETA ',filenam(1:3),' ',filenam(5:6),'   Power spectrum']);
                eval(['saveas(h2,''NO_THETA_ps_',filenam,'.fig'')']);
            case 2
                set(0,'CurrentFigure',h1);
                title(['SP THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
                eval(['saveas(h1,''SP_THETA_quo_',filenam,'.fig'')']);
                set(0,'CurrentFigure',h2);
                title(['SP THETA ',filenam(1:3),' ',filenam(5:6),'   Power spectrum']);
                eval(['saveas(h2,''SP_THETA_ps_',filenam,'.fig'')']);
            case 3
                set(0,'CurrentFigure',h1);
                title(['EV THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
                eval(['saveas(h1,''EV_THETA_quo_',filenam,'.fig'')']);
                set(0,'CurrentFigure',h2);
                title(['EV THETA ',filenam(1:3),' ',filenam(5:6),'   Power spectrum']);
                eval(['saveas(h2,''EV_THETA_ps_',filenam,'.fig'')']);
            end;
        end;
    end;
    waitbar(o/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
% fprf1 = ['k = ',num2str(k)];
% fprf2 = ['nmax = ',num2str(nmax)];
% fprf3 = ['m = ',num2str(m)];
% fprintf(fid,'%s\n%s\n%s',fprf1,fprf2,fprf3);
% fclose(fid);
close all;
cd(mm);

%   List of files calling RELFREQRUN: