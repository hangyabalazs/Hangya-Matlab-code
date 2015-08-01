function b_cont
%CONT   Plots the distribution of all theta segment length.
%   CONT loads the saved theta segments of all cells and calculates the distribution of
%   theta segment length. It saves the histogram in the result directory which you can
%   modify in the program code.
%
%   See also THETA and DISCONT.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where = [DATAPATH,'Wavelet\theta_segments\'];    %Here are the data files
files = dir(where);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH,'Wavelet\discont\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Loading THETA segments
allcont = [];
for o = 1:sf
    filenam = files(o).name;
    ff = [where filenam];
    load(ff)
    th = ThetaSegments;
    for k = 1:size(th,2)
        dis = th(2,k) - th(1,k);
        allcont(end+1) = dis;
    end
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator

% Plot histogram
[x,y] = hist(allcont,100);
B = bar(y,x);
eval(['saveas(B,''THETA_LENGTH_HIST.fig'')']);