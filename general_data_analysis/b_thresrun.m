function b_thresrun
%THRESRUN   Runs THRES on a sequence of files.
%   Automatic thresholding test. THRESRUN reads data from a 'where' directory and
%   saves the thresholds and the thresholded unit figures in a result directory.
%   You have to specify these directories in the program code.
%
%   See also THRES.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where = [DATAPATH,'DATA\analysenow2\'];    %Here are the data files
files = dir(where);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH,'DATA\thres\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator

for o = 1:sf
    fname = files(o).name;
    ffnm = [where fname];
    load(ffnm);
    if size(data,2) == 1
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
    
%Computing the input variables
    datinx1 = 1;
    datinx2 = 600000;
    b_imp(fname,where,data,datinx1,datinx2);
    
% Thresholding
    close all
    [T,H] = b_thres;
    kuszob = T;
   
% Saving
    flnm = fname(1:6);
    eval(['save(''THRES_',flnm,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'',''kuszob'')']);
    eval(['saveas(H,''THRES_',flnm,'_',num2str(datinx1),'_',num2str(datinx2),'.fig'')']);
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);