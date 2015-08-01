function b_thetaselectorrun2(inpdir,resdir)
%THETASELECTORRUN2   Runs THETASELECTOR3 and THETASELECTOR_BETA3 on a sequence of files
%   THETASELECTORRUN2 uses two directories: one for the data files and one for the results.
%   You are able to modify these directories through editing the program code, or you
%   can add them as input parameters: THETASELECTORRUN2(INPUT_DIRECTORY,RESULT_DIRECTORY).
%   Input directory can be specified using THETASELECTORRUN2(INPUT_DIRECTORY) syntax.
%
%   The function creates two subdirectories in the results directory: 'matrix' and 'jpg'.
%   It saves the output matrix (OM) in the 'matrix' subdirectory. OM is a 3-N matrix,
%   which contains THETASELECTOR3 output matrix (maximum locations and maximum values)
%   in the first and second row, and THETASELECTOR_BETA3 output matrix (theta-delta ratio)
%   in the third row. THETASELECTOR3 output figure (wavelet with maximum locations) gets
%   saved in 'jpg' subdirectory.
%
%   See also THETASELECTOR3, THETASELECTOR_BETA3 and THETASELECTORRUN.

% Input arguments check
error(nargchk(0,2,nargin));

% Define directories
global DATAPATH
if isempty(DATAPATH)
    clear global DATAPATH;
    b_startup
    global DATAPATH
end
global DATADIR
if isempty(DATADIR)
    clear global DATADIR;
    b_startup
    global DATADIR
end
if nargin < 2
    resdir = [DATAPATH, 'Wavelet2\'];  %Here are the results
end
if nargin < 1
    inpdir = [DATAPATH,'DATA\analysenow4\'];    %Here are the discriminated data files
end
cd(resdir)
create_subdir       %create subdirectories

% Import
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;

wb = waitbar(0,'Running THETASELECTORRUN2...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    fname = files(o).name;
    ffnm = [inpdir fname];
    eeg = [];
    load(ffnm);
    if isempty(eeg)
        eeg = data(:,1)';
    end
    
% Theta selection
    [pow,f,newstep] = powercall(eeg);
    [H,OM] = b_thetaselector3(pow,f,newstep);
    ratio = b_thetaselector_beta3(pow,f,newstep);
    
% Saving
    Out = zeros(3,size(pow,2));
    Out(1:2,:) = OM;
    Out(3,:) = ratio;
    
    pont = findstr(fname,'.');
    filenam = fname(1:pont(1)-1);
    fig = gcf;
    ax = findobj(fig,'Type','axes');
    axes(ax(2))
    fs = findstr(filenam,'_');
    fln1 = filenam(1:(fs(1)-1));
    fln2 = filenam((fs(1)+1):(fs(2)-1));
    title(['EEG ',fln1,' ',fln2]);
    cd jpg
    eval(['saveas(H,''THETA_SELECT_',fln1,'_',fln2,'.jpg'')']);
    cd ..
    cd matrix
    eval(['save(''THETA_SELECT_',fln1,'_',fln2,'.mat'',''Out'')']);
    cd ..
    clear pow
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all

% Save scalevector and 'newstep'
save scalevector f
save newstep newstep

% ----------------------------------------------------------------------------------
function create_subdir

% Create subdirectories
if ~b_isdir2('thetaselection')
    mkdir thetaselection
end
cd thetaselection
if ~b_isdir2('matrix')
    mkdir matrix
end
if ~b_isdir2('jpg')
    mkdir jpg
end

% ----------------------------------------------------------------------------------
function [pow,f,newstep] = powercall(wholeeeg)
%POWERCALL(WHOLEEEG)   Calls wavelet for EEG.
%   POWERCALL(WHOLEEEG) calls wavelet for EEG with fixed parameters (e.g. newstep for
%   downsampling is 25, which means 400 Hz downsampling using 10 kHz sampled data).
%
%   POWERCALL calls WAVELET_POW3.

% Create variables for wavelet
lenwe = length(wholeeeg);
newstep = 25;
resamp = 10000/newstep;
sst = wholeeeg(1:newstep:lenwe);

% Standardization
variance = std(sst)^2;
sst = (sst - mean(sst)) / sqrt(variance);

% Wavelet transformation
dt = 1 / resamp;    %resample on 400 Hz
time = [0:length(sst)-1] * dt + 0;
n = length(sst);
pad = 1;      
dj = 0.04;
s0 = 2 * dt;
j1 = ((1 / dj) * log2(n/2)) * 2;
j1 = ceil(j1);
j = (0:j1);
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;  
mother = 'Morlet';
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
[pow,period,scale,coi] = b_wavelet_pow3(sst,dt,pad,dj,s0,j1,mother,param,mis);