function b_megawaverun2
%MEGAWAVERUN2    Runs MEGAWAVE on a sequence of files.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   MEGAWAVERUN2 calculates and saves unit and EEG wavelet. It uses THRES automatic
%   thresholding function.
%
%   See also MEGAWAVERUN, MEGAWAVE and MEGAWAVEDISC.

% Input arguments check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
global DATADIR
where1 = [DATAPATH,'Data\analysenow4\'];    %Here are the data files
% where1 = DATADIR;
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH 'Entropy\megawave\']);  %Here are the results
if ~b_isdir2('thres')
    mkdir thres
end

% Import
wb = waitbar(0,'Running MEGAWAVERUN2...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
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
        datinx1 = 300000;
        datinx2 = 1200000;
        b_imp(fname,where1,data,datinx1,datinx2)
        
% Thresholding
        [kuszob,H] = b_thres;
        cd thres
        eval(['saveas(H,''THRES_',filenam,'_',num2str(datinx1),'_',num2str(datinx2),'.fig'')']);
        cd ..
        
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
end
close(wb);   %Close progress indicator
close all
cd(mm);