function b_ica_anal2
%ICA_ANAL2   Sorts out files if there is a decrease in the inter-first-spike interval variance.
%   This function uses three directories - one for the data files and two for the results. You 
%   have to specify these directories in the program code. The function sorts the ica-matrices 
%   into the two result directories. The sorting criterium is whether there is a decrease in the
%   inter-first-spike interval variance (examining only to 10 clusters).
%
%   See also ICA_BETA, ICA_BETA2, ICA_ANAL, ICA_ANAL3 and ICA_GUI.

% Input arguments check
error(nargchk(0,0,nargin));

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
respath1 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_decrease\'];  %Here are the results
respath2 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_nodecrease\'];  %Here are the results

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
    if ~files(o).isdir
        fname = files(o).name;
        ffnm = [where1 fname];
        load(ffnm)
        
% Computing the minimums and sorting the files
        dfsp = diff(FirstSpikeVar(3:11));
        ifd = isempty(find(dfsp<0));
        if ifd
            eval(['copyfile(''',ffnm,''',''',respath2,''');']);
        else
            eval(['copyfile(''',ffnm,''',''',respath1,''');']);
        end
    end
    waitbar(o/sf)   %Progress indicator
end

close(wb);   %Close progress indicator
cd(mm);