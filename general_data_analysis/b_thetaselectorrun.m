function long = b_thetaselectorrun(inpdir,resdir)
%THETASELECTORRUN   Runs THETASELECTOR3 and THETASELECTOR_BETA3 on a sequence of files
%   THETASELECTORRUN uses two directories: one for the data files and one for the results.
%   You are able to modify these directories through editing the program code, or you
%   can add them as input parameters: THETASELECTORRUN(INPUT_DIRECTORY,RESULT_DIRECTORY).
%   Input directory can be specified using THETASELECTORRUN(INPUT_DIRECTORY) syntax.
%
%   LONG = THETASELECTORRUN(INPDIR,RESDIR) returns file names with registrations longer
%   than 200 seconds in 'LONG' cell.
%
%   The function creates two subdirectories in the results directory: 'matrix' and 'jpg'.
%   It saves the output matrix (OM) in the 'matrix' subdirectory. OM is a 3-N matrix,
%   which contains THETASELECTOR3 output matrix (maximum locations and maximum values)
%   in the first and second row, and THETASELECTOR_BETA3 output matrix (theta-delta ratio)
%   in the third row. THETASELECTOR3 output figure (wavelet with maximum locations) gets
%   saved in 'jpg' subdirectory.
%
%   You have to use THETASELECTORRUN_LONG for registrations longer than 200 seconds.
%   THETASELECTORRUN saves a binary file called 'long.mat', which contains the identifiers
%   of the long registrations.
%
%   See also THETASELECTOR3, THETASELECTOR_BETA3 and THETASELECTORRUN_LONG.

% Input arguments check
error(nargchk(0,2,nargin));

% Define directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end
global DATADIR
if isempty(DATADIR)
    clear glogal DATADIR;
    b_startup
    global DATADIR
end
if nargin < 2
    resdir = [DATAPATH, 'Wavelet2\'];  %Here are the results
end
if nargin < 1
    inpdir = [DATAPATH,'DATA\analysenow4\'];    %Here are the data files
end
cd(resdir)
create_subdir       %create subdirectories

% Dealing with long files - clear global 'LONG'
clear glogal LONG;

% Import
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;

wb = waitbar(0,'Running THETASELECTORRUN...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf,
    fname = files(o).name;
    ffnm = [inpdir fname];
    data = b_load_data(ffnm);
    meret = size(data,1);
    
% Computing the input variables
    datinx1 = 1; %first point of the interval
    datinx2 = size(data,1); %last point of the interval
    if datinx2 > 2000000
        global LONG
        LONG{end+1} = fname;
        waitbar(o/sf)   %Progress indicator
        continue
    end
    b_imp(fname,inpdir,data,datinx1,datinx2)
    
% Free memory
    save temp o sf wb files inpdir
    clear
    clear function
    
% Theta selection
    [wave,f,newstep] = b_waveletcall_for_applications;
    
    sw2 = size(wave,2);
    pieceno = 5;
    segm = fix(sw2/pieceno);
    power = [];
    while ~isempty(wave)
        index1 = 1;
        index2 = min(segm,size(wave,2));
        wavefrag = wave(:,index1:index2);
        powerfrag = (abs(wavefrag)) .^ 2;
        clear wavefrag
        wave(:,index1:index2) = [];
        power = [power powerfrag];
    end
    clear wave
    
    [H,OM] = b_thetaselector3(power,f,newstep);
    ratio = b_thetaselector_beta3(power,f,newstep);
    
% Saving
    Out = zeros(3,sw2);
    Out(1:2,:) = OM;
    Out(3,:) = ratio;
    
    global IN
    fname = IN{3};
    pont = findstr(fname,'.');
    filenam = fname(1:pont(1)-1);
    fig = gcf;
    ax = findobj(fig,'Type','axes');
    axes(ax(2))
    title(['EEG ',filenam(1:3),' ',filenam(5:6)]);
    cd jpg
    eval(['saveas(H,''THETA_SELECT_',filenam(1:6),'.jpg'')']);
    cd ..
    cd matrix
    eval(['save(''THETA_SELECT_',filenam(1:6),'.mat'',''Out'')']);
    cd ..
    load temp
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
delete temp.mat
close all

% Save 'long'
global LONG
long = LONG;
save long long
clear global LONG

% Save scalevector and 'newstep'
save scalevector f
save newstep newstep

% ----------------------------------------------------------------------------------
function create_subdir

% Create subdirectories
if ~b_isdir2(thetaselection)
    mkdir(thetaselection)
end
cd thetaselection
if ~b_isdir2('matrix')
    mkdir matrix
end
if ~b_isdir2('jpg')
    mkdir jpg
end