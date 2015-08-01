% --------------------------------------------------------------------------------------
% Modify the function name and the help!
% --------------------------------------------------------------------------------------
function b_empty
%EMPTY    Scheme for programs running on sequence of files.
%   EMPTY is an uncomplete program which you can use for batch processing on hippocampus
%   data files. You have to write the applied program in the code. See the comments for
%   help.
%
%   See also EMPTY_NEW, EMPTY_2004 and EMPTY_GUI2004..

% Input arguments check
% --------------------------------------------------------------------------------------
% Replace first zero with the minimum and second zero with the maxaimum number of input 
% arguments, if there is any!
% --------------------------------------------------------------------------------------
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
% --------------------------------------------------------------------------------------
% This function uses two directories - one for the data files and one for the results.
% You have to specify these directories in the following lines.
% --------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------
% Replace *** with the name of the directory containing the data files!
% --------------------------------------------------------------------------------------
where1 = [DATAPATH,'***'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
% --------------------------------------------------------------------------------------
% Replace *** with the name of the directory in which you want to put the results!
% --------------------------------------------------------------------------------------
cd([DATAPATH,'***']);  %Here are the results
% --------------------------------------------------------------------------------------
% Comment the following 18 lines out, if you do not want to create text file!
% --------------------------------------------------------------------------------------
npt = input...
% --------------------------------------------------------------------------------------
% Replace *** with the name of the text file you want to create!
% --------------------------------------------------------------------------------------
    ('Discard existing content or append data while writing ***? /discard:  ENTER, append: a/','s');
if isempty(npt)
% --------------------------------------------------------------------------------------
% Replace *** with the name of the text file you want to create!
% --------------------------------------------------------------------------------------
    fid = fopen('***','w');
elseif npt == 'a'
% --------------------------------------------------------------------------------------
% Replace *** with the name of the text file you want to create!
% --------------------------------------------------------------------------------------
    fid = fopen('***','a');
else
    error('Unexpected answer for input.')
end

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = b_load_data(ffnm);
    
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
            b_imp(fname,where1,data,datinx1,datinx2)
            
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
% --------------------------------------------------------------------------------------
% If you do not create a text file, comment out the next line!
% If you create one, you have to insert the fprint lines in the code!
% --------------------------------------------------------------------------------------
fclose(fid);

close all
cd(mm);