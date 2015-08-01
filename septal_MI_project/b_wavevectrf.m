function b_wavevectrf
%WAVEVEC_TRF    Rename files.

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'ICA\ica_wave3\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
pathname = [DATAPATH,'analysenow1\'];    %Here are the datinx 
mm = pwd;
cd([DATAPATH,'ICA\ica_wave4\']);  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    fname = files(o).name;
    ffnm = [where1 fname];
    load(ffnm);
    i_first = fname(12:17);
    i_second = fname(19:24);
    cl = fname(26:31);
    str = ['WAVEVECTOR_',cl,'_',i_first,'_',i_second];
    eval(['save ' str ' ' 'Wavevec']);
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator