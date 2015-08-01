function b_thresrun2_3ch
%THRESRUN2_3CH   Runs automatic thresholder on a sequence of 3-channel data files.
%   THRESRUN2_3CH uses two directories: one for the input files and another for the results;
%   you are able to modify these directories through editing the program code.
%
%   THRESRUN2_3CH saves threshold values in the results' directory. After
%   running THRESRUN2_3CH you are able to create new type data files 
%   (containing eeg and vdisc) using THRES_GUI.
%
%   THRESRUN2_3CH is a version of THRESRUN2 implemented for 3-channel data
%   files.
%
%   See also THRESRUN2, THRES5, DISC2, DISCRUN and THRES_GUI.

% Input arguments check
error(nargchk(0,0,nargin));

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR '3ch\temp\'];    %Here are the data files
resdir = [DATAPATH,'MI\Thres_hc\'];            %Here are the results
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;
cd(resdir);

% Import
wb = waitbar(0,'Running THRESRUN2 3CH...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    fname = files(o).name;
    filenam = fname(1:6);
    ffnm = [inpdir fname];
    data = b_load_data(ffnm);
    
% Computing the input variables
    datinx1 = 1;
    datinx2 = length(data);
    for i = 1:2
        if isequal(i,1)
            HCimp(fname,inpdir,data,datinx1,datinx2);   % i = 1: HC unit
        else
            MSimp(fname,inpdir,data,datinx1,datinx2);   % i = 2: MS unit
        end
        
% Thresholding
        close all
        [T,seglen] = b_thres5;
        
% Discrimination
        b_disc2(T,seglen);
        global DISC
        id = DISC{1};
        output = DISC{2};
        vdisc = DISC{3};

% Saving
        if isequal(i,1)
            str = ['THRES_HCunit_',fname];
        else
            str = ['THRES_MSunit_',fname];
        end
        save(str,'T','seglen','vdisc');
    end
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
close all
cd(mm);

% ------------------------------------------------------------------------
function HCimp(fname,where,data,datinx1,datinx2);

% Creating the import variables
unit = data(datinx1:datinx2,1);
unit = unit';
eeg = data(datinx1:datinx2,3);
eeg = eeg';
dt = 0.0001;
time = [0:length(unit)-1] * dt; 
xlimit = [min(time),max(time)];
meret = size(data,1);
mintafr = 1 / dt;

% Create global IN
global IN
IN = cell(1,12);
IN{1} = data;
IN{2} = eeg;
IN{3} = fname;
IN{4} = where;   %path name
IN{5} = datinx1;
IN{6} = datinx2;
IN{7} = time;
IN{8} = unit;
IN{9} = dt;
IN{10} = meret;
IN{11} = mintafr;
IN{12} = xlimit;

% Assigning the variables in the caller workspace
ws = 'caller';
assignin(ws,'eeg',IN{2})
assignin(ws,'time',IN{7})
assignin(ws,'unit',IN{8})
assignin(ws,'dt',IN{9})
assignin(ws,'meret',IN{10})
assignin(ws,'mintafr',IN{11})
assignin(ws,'xlimit',IN{12})

% ------------------------------------------------------------------------
function MSimp(fname,where,data,datinx1,datinx2);

% Creating the import variables
unit = data(datinx1:datinx2,2);
unit = unit';
eeg = data(datinx1:datinx2,3);
eeg = eeg';
dt = 0.0001;
time = [0:length(unit)-1] * dt; 
xlimit = [min(time),max(time)];
meret = size(data,1);
mintafr = 1 / dt;

% Create global IN
global IN
IN = cell(1,12);
IN{1} = data;
IN{2} = eeg;
IN{3} = fname;
IN{4} = where;   %path name
IN{5} = datinx1;
IN{6} = datinx2;
IN{7} = time;
IN{8} = unit;
IN{9} = dt;
IN{10} = meret;
IN{11} = mintafr;
IN{12} = xlimit;

% Assigning the variables in the caller workspace
ws = 'caller';
assignin(ws,'eeg',IN{2})
assignin(ws,'time',IN{7})
assignin(ws,'unit',IN{8})
assignin(ws,'dt',IN{9})
assignin(ws,'meret',IN{10})
assignin(ws,'mintafr',IN{11})
assignin(ws,'xlimit',IN{12})