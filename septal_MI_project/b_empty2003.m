% --------------------------------------------------------------------------------------
% Modify the function name and the help!
% --------------------------------------------------------------------------------------
function b_empty2003(DataDir,ResultDir,IntervalSet,TextOutput,EvSpNo,InputArg)
%EMPTY2003    Scheme for programs running on sequence of files.
%   EMPTY2003 is an uncomplete program which you can use for batch processing on hippocampus
%   data files. You have to write the applied program in the code. See the comments for
%   help.
%
%   Different from EMPTY, you have to specify some parameters via EMPTY_GUI2003 Graphical
%   User Interface. Start EMPTY_GUI2003 to run!
%
%   See also EMPTY, EMPTY_NEW, EMPTY2004, EMPTY_GUI2003 and EMPTY_GUI2004.

% Input arguments
if ~isempty(InputArg)
    sia = findstr(InputArg,' ');
    sia = [0 sia length(InputArg)+1];
    arg = cell(1,length(sia));
    for i = 2:length(sia)
        inp = InputArg(sia(i-1)+1:sia(i)-1);
        if isnumeric(inp)
            arg{i} = inp;
        else
            arg{i} = evalin('base','eval(inp)');
        end
    end
end

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end

where1 = DataDir;    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
cd([DATAPATH,ResultDir]);  %Here are the results
if ~isempty(TextOutput)
    npt = input(['Discard existing content or append data while writing' ...
            TextOutput '? /discard:  ENTER, append: a/'],'s');
    if isempty(npt)
        str = ['fid = fopen(' TextOutput ',''w'');'];
        eval(str);
    elseif npt == 'a'
        str = ['fid = fopen(' TextOutput ',''a'');'];
        eval(str);
    else
        error('Unexpected answer for input.')
    end 
end

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
    end
    if size(data,2)==1,
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
    
%Searching for the matching datinx file
    if IntervalSet == 'auto'
        cl = fname(1:6);
        filename = [cl,'.txt'];
        pont = findstr(filename,'.');
        filenam = filename(1:pont(1)-1);
        fn = fullfile(pathname,filename);
        [s,v] = textread(fn,'%s%f');
        ww = 0;
    
        newp = [];
        if EvSpNo(1)
            newp = 1;
        end
        if EvSpNo(2)
            newp = [newp 2];
        end
        if EvSpNo(1)
            newp = [newp 3];
        end
    else
        newp = 0;
    end
    
            
%Computing the input variables
    for p = 1:length(newp)
        if newp(1) ~= 0
            datinx1 = v(2*newp(p)-1); %first point of the interval
            datinx2 = v(2*newp(p)); %last point of the interval
        else
            datinx1 = input('First point of the interval:');
            datinx2 = input('Last point of the interval:');
        end
   
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
            
% --------------------------------------------------------------------------------------
% This is the place of the computation!
% --------------------------------------------------------------------------------------
            
% --------------------------------------------------------------------------------------
% You have to modify or comment out the saving part, otherwise you get an error message!
% --------------------------------------------------------------------------------------
%Saving
            pont = findstr(filename,'.');
            filenam = filename(1:pont(1)-1);
            switch p
            case 1
                title(['NO THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
                eval(['saveas(h1,''NO_THETA_quo_',filenam,'.fig'')']);
            case 2
                title(['SP THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
                eval(['saveas(h1,''SP_THETA_quo_',filenam,'.fig'')']);
            case 3
                title(['EV THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
                eval(['saveas(h1,''EV_THETA_quo_',filenam,'.fig'')']);
            end
        end
    end
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
if ~isempty(TextOutput)
    fclose(fid);
end

close all
cd(mm);