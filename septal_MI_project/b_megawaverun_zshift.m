function b_megawaverun_zshift
%MEGAWAVERUN_ZSHIFT    Runs MEGAWAVE_ZSHIFT on a sequence of files.
%   This function uses two directories - one for the data files and one for
%   the results. You have to specify these directories in the program code.
%
%   MEGAWAVERUN_ZSHIFT calculates and saves eeg and shifted unit wavelet. 
%   Z-shift is chosen as 'delay' (see ZSHIFT for details). It uses 
%   discriminated input files.
%
%   See also MEGAWAVERUN2, MEGAWAVERUN_ENTRCTRL2, ZSHIFT and MEGAWAVE_ZSHIFT.

% Input arguments check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH
global DATADIR
where1 = [DATAPATH,'Data\analysenow4\'];    %Here are the data files
% where1 = DATADIR;
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH 'Entropy7\megawave\']);  %Here are the results

% Import
wb = waitbar(0,'Running MEGAWAVE ZSHIFT...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

sf = min(11,sf);
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
    [wavea,wavez] = b_megawave_zshift(eeg,vdisc,filenam);
    
% Saving
    filnm = ['MEGAWAVE_RESULT_' filenam '.mat'];
    save(filnm,'wavea','wavez')
    clear wavea wavez
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);