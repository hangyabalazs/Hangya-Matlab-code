function b_discrun
%DISCRUN   Runs automatic thresholder and discriminator on a sequence of files.
%   DISCRUN uses two directories: one for the input files and another for the results;
%   you are able to modify these directories through editing the program code.
%
%   DISCRUN saves threshold values and 'vdisc' in the results' directory. After
%   running DISCRUN you are able to create new type data files (containing eeg 
%   and vdisc) using THRES_GUI.
%
%   See also THRES4, DISC2, THRESRUN2 and THRES_GUI.

% Input arguments check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
global DATADIR
% inpdir = [DATAPATH,'DATA\analysenow2\'];    %Here are the data files
inpdir = [DATADIR 'PViktor\'];
resdir = [DATAPATH 'PViktor\thres\'];            %Here are the results
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;
cd(resdir);

% Import
wb = waitbar(0,'Running DISCRUN...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    fname = files(o).name;
    filenam = fname(1:6);
    ffnm = [inpdir fname];
    data = b_load_data(ffnm);
    
% Computing the input variables
    datinx1 = 1;
    datinx2 = length(data);
    b_imp(fname,inpdir,data,datinx1,datinx2);
    
% Thresholding
    close all
    [T,seglen] = b_thres4;
    
% Discrimination
    b_disc2(T,seglen);
    global DISC
    id = DISC{1};
    output = DISC{2};
    vdisc = DISC{3};
    
    str = ['THRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'];
    eval(['save ' str ' T seglen vdisc']);
    
% Saving
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);