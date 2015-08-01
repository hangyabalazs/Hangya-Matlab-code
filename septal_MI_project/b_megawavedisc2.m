function b_megawavedisc2
%MEGAWAVEDISC2   Saves thresholds for unit discrimination.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   The differences between MEGAWAVEDISC and MEGAWAVEDISC2: 'where1', 'datinx1', 'datinx2'.
%
%   See also DISC, MEGAWAVEDISC and MEGAWAVERUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'Data\analysenow2\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH,'Data\megawavedisc3\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    data = load(ffnm);
    meret=size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end;
    if size(data,2)==1,
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end;
    
% Computing the input variables
    
    datinx1 = 320000; %first point of the interval
    datinx2 = 920000; %last point of the interval
    
    dt = 0.0001;
    mintafr = 1 / dt;
    unit = data(datinx1:datinx2,2);
    unit = unit';
    eeg = data(datinx1:datinx2,1);
    eeg = eeg';
    time = [0:length(unit)-1] * dt; 
    xlimit = [min(time),max(time)];
    
% Finding the threshold                
    s = figure;
    plot(unit,'m');
    set(gca,'xlim',[0 length(unit)])
    title('Give the threshold! /Command window/')
    kuszob = input('Give the threshold! ');
    
% Saving
    eval(['save(''MEGAWAVE_',filenam,'.mat'',''kuszob'')']);
    
    waitbar(o/sf)   %Progress indicator
    close all
end
close(wb);   %Close progress indicator

close all