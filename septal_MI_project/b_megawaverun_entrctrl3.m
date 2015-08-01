function b_megawaverun_entrctrl3
%MEGAWAVERUN_ENTRCTRL3    Runs MEGAWAVE_ENTRCLTRL3 on a sequence of files.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   MEGAWAVERUN_ENTRCTRL3 calculates and saves delayed unit wavelet (with 10 seconds delay).
%   It uses discriminated input files.
%
%   See also MEGAWAVERUN2, MEGAWAVERUN_ENTRCTRL and MEGAWAVE_ENTRCTRL3.

% Input arguments check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH
global DATADIR
where1 = [DATAPATH,'Data\analysenow3\'];    %Here are the data files
% where1 = DATADIR;
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH 'Entropy\control\megawave2\']);  %Here are the results

% Import
wb = waitbar(0,'Running MEGAWAVE ENTRCTRL3...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

sf = min(10,sf);
for o = 1:sf
    fname = files(o).name;
    ffnm = [where1 fname];
    filenam = fname(1:6);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    load(ffnm);
    
% Computing the input variables
    datinx1 = 300000;
    datinx2 = 1200000;
    eeg = eeg(datinx1:datinx2);
    vdisc = vdisc(find(vdisc>=datinx1&vdisc<=datinx2)) - datinx1;
    
% Wavelet transformation
    [wavef] = b_megawave_entrctrl3(eeg,vdisc);
    
% Saving
    filnm = ['MEGAWAVE_RESULT_' filenam '.mat'];
    save(filnm,'wavef')
    clear wavef
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);