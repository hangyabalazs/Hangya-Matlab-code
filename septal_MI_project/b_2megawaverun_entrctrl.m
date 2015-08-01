function b_2megawaverun_entrctrl
%2MEGAWAVERUN_ENTRCTRL    Runs 2MEGAWAVE_ENTRCLTRL on a sequence of files.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   2MEGAWAVERUN_ENTRCTRL calculates and saves eeg, unit, random sample, random generated
%   unit and delayed unit wavelet. It uses discriminated input files.
%
%   See also MEGAWAVERUN2 and 2MEGAWAVE_ENTRCTRL.

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
cd([DATAPATH 'Entropy2\control\megawave\']);  %Here are the results

% Import
wb = waitbar(0,'Running 2MEGAWAVE ENTRCTRL...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

sf = min(31,sf);
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
    [wavea,waveb,wavec,waved,wavee,wavef] = b_2megawave_entrctrl(eeg,vdisc);
    
% Saving
    filnm = ['MEGAWAVE_RESULT_' filenam '.mat'];
    save(filnm,'wavea','waveb','wavec','waved','wavee','wavef')
    clear wavea waveb wavec waved wavee wavef
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);