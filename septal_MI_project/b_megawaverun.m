function b_megawaverun
%MEGAWAVERUN    Runs MEGAWAVE on a sequence of files.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   MEGAWAVERUN calculates and saves unit and EEG wavelet. It loads previously saved
%   threshold values, so you should be careful with modifying the segment limits!
%
%   See also MEGAWAVERUN2, MEGAWAVE_FOR_MEGAWAVERUN and MEGAWAVEDISC.

% Input arguments check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'Data\analysenow3\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH 'Entropy\megawave']);  %Here are the results

% Import
wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    data = b_load_data(ffnm);
    
% Computing the input variables
%     segmlen = 80000;    %we analyse the first and the second minute of the data after each other 
%     for p = 4:12
%         datinx1 = (p-1) * segmlen + 1; %first point of the interval
%         datinx2 = p * segmlen; %last point of the interval
        datinx1 = 500000;
        datinx2 = 700000;
        b_imp(fname,where1,data,datinx1,datinx2)
        
% Loading the saved threshold value
        directory = [DATAPATH,'Data\megawavedisc3\'];
        str = [directory,'MEGAWAVE_','_',filenam,'.mat'];
        load(str)
        
% Discrimination                
        b_disc(kuszob);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};
        kuszob = DISC{4};
        instfrek = DISC{5};
        isi = DISC{6};
        
% Wavelet transformation
        [wavea,waveb] = b_megawave;
% Saving
%         eval(['save(''MEGAWAVE_RESULT_','_',filenam,'.mat'',''wavea'',''waveb'')']);
          filnm = ['MEGAWAVE_RESULT_' filenam '.mat'];
          save(filnm,'wavea','waveb')
          clear wavea waveb
%     end
    waitbar(o/sf)   %Progress indicator
    close all
end
close(wb);   %Close progress indicator
% fclose(fid);

close all
cd(mm);