function b_thetaicatrf
%THETAICATRF    Rename files.

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
where2 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly2\'];  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf,
    if files(o).isdir
        continue
    end
    fname = files(o).name;
    i_first = fname(11:16);
    i_second = fname(18:23);
    cl = fname(25:30);
    fname2 = ['THETA_ICA_',cl,'_',i_first,'_',i_second,'.mat'];
    source = [where1 fname];
    dest = [where2 fname2];
    eval(['copyfile ' source  ' ' dest]);
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator