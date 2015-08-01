function b_perioddrun
%PERIODRUN   Runs LOMB_PERIOD on a sequence of files.
%   This function uses two directories - one for the data files and one for the results.
%   You have to specify these directories in the program code.
%
%   PERIODRUN creates .mat files for LOMBPER_FOR_ICA which is used by ICA_GUI2B. Lomb-
%   Scargle periodogram is calculated for "burst all spikes": sampling points are the
%   locations of intraburst spikes and sampled data is the Blackmann-Harris convolved
%   unit.
%
%   See also ICA_GUI2B, LOMBPER_FOR_ICA, LOMB_PERIOD and LOMBPER.

% Input arguments check
error(nargchk(0,0,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'ICA\ica_beta_thetaonly2\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
% cd([DATAPATH,'ICA\lomb\']);  %Here are the results
% npt = input...
%     ('Discard existing content or append data while writing ***? /discard:  ENTER, append: a/','s');
% if isempty(npt)
%     fid = fopen('***','w');
% elseif npt == 'a'
%     fid = fopen('***','a');
% else
%     error('Unexpected answer for input.')
% end

% wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    load(ffnm);
    i_first = str2num(fname(11:16));
    i_second = str2num(fname(18:23));
    
% LOMB periodogram
    Vdisc_new = Vdisc - (i_first - 320000);
    d = length(IntraBurstIvVar);
    for dec = 3:d
        bas = [];
        for bno = 1:size(Burst{dec},2)
            b = Vdisc_new(Burst{dec}(1,bno):Burst{dec}(2,bno));
            bas = [bas b];
        end
        lenu = length(Time);
        vdisc_for_lomb = Vdisc_new;
        [pxx,pyy,jmax,prob,z,effm,period_handle] = b_lomb_period(vdisc_for_lomb,bas,lenu);
        
%Saving
        pont = findstr(fname,'.');
        filenam = fname(1:pont(1)-1);
        flnm = filenam(25:30);
%         eval(['save(''LOMB_',flnm,'_',num2str(dec),'_',num2str(i_first),'_',num2str(i_second),'.mat'',''pxx'',''pyy'',''z'')']);
    end
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
% fclose(fid);
close all
cd(mm);