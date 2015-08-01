function b_discont
%DISCONT   Plots the distribution of all theta discontinuities.
%   DISCONT loads the saved theta segments of all cells and calculates the distribution of
%   inter-theta segment length. It saves the histogram in the result directory which you
%   can modify in the program code.
%
%   See also THETA and CONT.

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
alldis = [];
for o = 1:sf
    filenam = files(o).name;
    ff = [where filenam];
    load(ff)
    th = ThetaSegments;
    for k = 1:size(th,2)-1
        dis = th(1,k+1) - th(2,k);
        alldis(end+1) = dis;
    end
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator

% Plot histogram
[x,y] = hist(alldis,100);
B = bar(y,x);
eval(['saveas(B,''DISCONTINUITY_HIST.fig'')']);