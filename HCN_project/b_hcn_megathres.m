function b_hcn_megathres
%HCN_MEGATHRES   Saves thresholds for unit discrimination.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   HCN_MEGATHRES saves the "long file" names (i.e. registration is longer than 200 sec.)
%   in the results' directory. You have to threshold them with HCN_MEGATHRES_LONG.
%
%   See also DISC, DISC2, MEGAWAVEDISC, MEGAWAVERUN and HCN_MEGATHRES_LONG.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = ['f:\raw_data\hcn\temp\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
    end
end
files2 = files2(2:end);
sf = length(files2);
mm = pwd;
cd([DATAPATH,'HCN\Data\thresholds\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
    fname = files2(o).name;
    ffnm = [where1 fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    data = load(ffnm);
    meret = size(data,1);
    if isstruct(data),
        field = fieldnames(data);
        s = struct('type','.','subs',field);
        data = subsref(data,s);
    end
    if size(data,2)==1
        data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
    end
    
% Computing the input variables
    datinx1 = 1;     %first point of the interval
    datinx2 = size(data,1);      %last point of the interval
    if datinx2 > 2000000
        global LONG
        LONG{end+1} = fname(1:6);
        waitbar(o/sf)   %Progress indicator
        continue
    end
    
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
    eval(['save(''MEGATHRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.mat'',''kuszob'')']);
    
    waitbar(o/sf)   %Progress indicator
    close all
end
if exist(LONG)
    save long LONG
end
close(wb);   %Close progress indicator
close all