function b_sdffft
%SDFFFT Calculates fft of unit spike density function on a sequence of files.
%   This program is to find objectiv criteria of theta frequency bursting (working 
%   on burst_cls_auto).
%
%   The function uses three directories - one for the data files, one for the datinx
%   (first and last points of "no theta", "spontanous theta" and "evoked theta" intervals,
%   thresholds for the discrimination) and one for the results. You have to specify these
%   directories in the program code. The function saves the created plot in the result
%   directory. It also saves the relative and the maximum theta power in a text file.
%
%   See also FFT.

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
cd([DATAPATH,'sdffft_figs_proba\']);  %Here are the results
npt = input...
    ('Discard existing content or append data while writing sdffft.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('sdffft.txt','w');
elseif npt == 'a',
    fid = fopen('sdffft.txt','a');
else
    error('Unexpected answer for input.')
end;

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for i = 1:sf,
    fname = files(i).name;
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
   
%Searching for the matching datinx file
    cl = fname(1:6);
    filename = [cl,'.txt'];
    pont = findstr(filename,'.');
    filenam = filename(1:pont(1)-1);
    fn = fullfile(pathname,filename);
    [s,v] = textread(fn,'%s%f');
    ww = 0;

%Computing the input variables
    for i = 1:3,
        datinx1 = v(2*i-1); %first point of the interval
        datinx2 = v(2*i); %last point of the interval
   
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
    

%Computing Spike Density Function
            ipunit = zeros(1,length(unit)); %SDF
            ipunit(vdisc) = 1;
            wbh = blackmanharris(1000);
            wipunit = conv(ipunit,wbh);
            wipunit = wipunit(1:length(unit));
            [punit funit] = periodogram(wipunit,blackmanharris(length(wipunit)),65536,10000);
            relthetapower = (sum(punit(21:41))) / (sum(punit(4:end)));    %relative theta power
            maxamp = max(punit(21:41));    %maximum  theta power
            
%Plotting and saving
            pont = findstr(filename,'.');
            filenam = filename(1:pont(1)-1);
            switch i
            case 1
                subplot(2,1,1);
                h1 = plot(funit(1:100),punit(1:100));
                axis([0 16 0 2]);
                text(8,1.5,['{\itRelative theta power: }',num2str(relthetapower)]);
                text(8,1,['{\itMaximum theta power: }',num2str(maxamp)]);
                subplot(2,1,2);
                plot(time,unit,'m');
                title(['NO THETA ',filenam(1:3),' ',filenam(5:6)]);
                eval(['saveas(h1,''NO_THETA_',filenam,'.fig'')']);
                ffilenam = [filenam,' NO THETA'];
                fprf = [ffilenam,' ', num2str(relthetapower),' ', num2str(maxamp)];
                fprintf(fid,'%s %s %s\n',fprf);
                fprintf(fid,'\n',fprf);
            case 2
                subplot(2,1,1);
                h2 = plot(funit(1:100),punit(1:100));
                axis([0 16 0 2]);
                text(8,1.5,['{\itRelative theta power: }',num2str(relthetapower)]);
                text(8,1,['{\itMaximum theta power: }',num2str(maxamp)]);
                subplot(2,1,2);
                plot(time,unit,'m');
                title(['SP THETA ',filenam(1:3),' ',filenam(5:6)]);
                eval(['saveas(h2,''SP_THETA_',filenam,'.fig'')']);
                ffilenam = [filenam,' SP THETA'];
                fprf = [ffilenam,' ', num2str(relthetapower),' ', num2str(maxamp)];
                fprintf(fid,'%s %s %s\n',fprf);
                fprintf(fid,'\n',fprf);
            case 3
                subplot(2,1,1);
                h3 = plot(funit(1:100),punit(1:100));
                axis([0 16 0 2]);
                text(8,1.5,['{\itRelative theta power: }',num2str(relthetapower)]);
                text(8,1,['{\itMaximum theta power: }',num2str(maxamp)]);
                subplot(2,1,2);
                plot(time,unit,'m');
                title(['EV THETA ',filenam(1:3),' ',filenam(5:6)]);
                eval(['saveas(h3,''EV_THETA_',filenam,'.fig'')']);
                ffilenam = [filenam,' EV THETA'];
                fprf = [ffilenam,' ', num2str(relthetapower),' ', num2str(maxamp)];
                fprintf(fid,'%s %s %s\n',fprf);
                fprintf(fid,'\n',fprf);
            end;
        end;
    end;
    waitbar(i/sf)   %Progress indicator
end;
close(wb);   %Close progress indicator
fclose(fid);
cd(mm);

%   List of files calling SDFFFT: