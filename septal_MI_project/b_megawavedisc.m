function b_megawavedisc
%MEGAWAVEDISC   Saves thresholds for unit discrimination.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   See also DISC, MEGAWAVEDISC2 and MEGAWAVERUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'Data\analysenow2b\'];    %Here are the data files
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
%     segmlen = 80000;    %we analyse the first and the second minute of the data after each other 
%     for p = 4:12
%         datinx1 = (p-1) * segmlen + 1; %first point of the interval
%         datinx2 = p * segmlen; %last point of the interval
        datinx1 = 500000;
        datinx2 = 700000;
        
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
        
%     end
    waitbar(o/sf)   %Progress indicator
    close all
end
close(wb);   %Close progress indicator
% fclose(fid);

close all
% cd(mm);