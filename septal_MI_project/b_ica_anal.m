function b_ica_anal(n)  %n=0.1
%ICA_ANAL   Sorts files by minimum of inter-first-spike interval variance.
%   This function uses three directories - one for the data files and two for the results.
%   You have to specify these directories in the program code. The function sorts the 
%   ica-matrices into the two result directories by the minimum of the 3rd to 8th values
%   of InterFirstSpikeVar matrix, spliting at the value of the input variable n.
%
%   See also ICA_BETA, ICA_BETA2, ICA_ANAL2, ICA_ANAL3 and ICA_GUI.

% Input arguments check
error(nargchk(1,1,nargin));

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
respath1 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_under\'];  %Here are the results
respath2 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_over\'];  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
    if ~files(o).isdir
        fname = files(o).name;
        ffnm = [where1 fname];
        load(ffnm)
        
% Computing the minimums and sorting the files
        mfsp = min(FirstSpikeVar(3:8));
        if mfsp > n
            eval(['copyfile(''',ffnm,''',''',respath2,''');']);
        else
            eval(['copyfile(''',ffnm,''',''',respath1,''');']);
        end
    end
    waitbar(o/sf)   %Progress indicator
end

close(wb);   %Close progress indicator
cd(mm);