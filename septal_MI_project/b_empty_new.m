function b_empty_new
%EMPTY_NEW    Scheme for programs running on sequence of files.
%   EMPTY_NEW is an uncomplete program which you can use for batch processing on hippocampus
%   data files. You have to write the applied program in the code. See the comments for
%   help.
%
%   See also EMPTY, EMPTY_2004 and EMPTY_GUI2004.

% --------------------------------------------------------------------------------------
% NAME & HELP:
%       Modify the function name and the help!
%
% INPUT ARGUMENTS CHECK
%       Replace first zero with the minimum and second zero with the maxaimum number of
%       input arguments, if there is any!
%
% DIRECTORIES:
%       This function uses two directories - one for the data files and one for the 
%       results. You have to specify these directories in the following lines.
%       Replace *** with the name of the directory!
%
% TEXT FILE:
%       In case you want to create a text file set istext to 1, otherwise istext is 0.
%       Replace *** with the name of the text file!
%
% MAIN and SAVE
%       You have to insert the call of the main function and modify the saving part at 
%       the end of the code.
% --------------------------------------------------------------------------------------

% Input arguments check
error(nargchk(0,0,nargin));

% Directories
where1 = [DATAPATH,'***'];    %Here are the data files
cd([DATAPATH,'***']);  %Here are the results

% Text File
istext = 0;
if istext
    txtfile = '***';
end

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;

if istext
    npt = input...
        (['Discard existing content or append data while writing' txtfile '? /discard:  ENTER, append: a/','s']);
    if isempty(npt)
        fid = fopen(txtfile,'w');
    elseif npt == 'a'
        fid = fopen(txtfile,'a');
    else
        error('Unexpected answer for input.')
    end
end

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    data = b_load_data(ffnm);
    
% Computing the input variables
    datinx1 = 1; %first point of the interval
    datinx2 = size(data,1); %last point of the interval
    b_imp(fname,where,data,datinx1,datinx2);
             
% --------------------------------------------------------------------------------------
% This is the place of the computation!
% --------------------------------------------------------------------------------------
            
% --------------------------------------------------------------------------------------
% You have to modify or comment out the saving part, otherwise you get an error message!
% --------------------------------------------------------------------------------------
%Saving
    pont = findstr(filename,'.');
    filenam = filename(1:pont(1)-1);
    title(['NO THETA ',filenam(1:3),' ',filenam(5:6),'   Quociens']);
    eval(['save(''THRES_',flnm,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'',''kuszob'')']);
    eval(['saveas(H,''THRES_',flnm,'_',num2str(datinx1),'_',num2str(datinx2),'.fig'')']);
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
if istext
    fclose(fid);
end

close all
cd(mm);