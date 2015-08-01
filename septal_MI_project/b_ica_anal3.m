function v = b_ica_anal3(n)  %n=0.15
%ICA_ANAL3   Sorts files by minimum of intraburst interval variance.
%   This function uses three directories - one for the data files and two for the results.
%   You have to specify these directories in the program code. The function sorts the 
%   ica-matrices into the two result directories by the 3rd value of IntraBurstIvVar 
%   matrix, spliting at the value of the input variable n. ICA_ANAL3 consideres only
%   the files in ica_beta_under directoy (previously sorted by ICA_ANAL).
%
%   See also ICA_BETA, ICA_BETA2, ICA_ANAL, ICA_ANAL2 and ICA_GUI.

% Input arguments check
error(nargchk(1,1,nargin));

% Import
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end;
where1 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_under\'];    %Here are the data files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
respath1 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_intraunder\'];  %Here are the results
respath2 = [DATAPATH,'ICA\ica_gui2b\ica_beta_thetaonly\ica_beta_intraover\'];  %Here are the results
npt = input...
    ('Discard existing content or append data while writing anal_3.txt? /discard:  ENTER, append: a/','s');
if isempty(npt),
    fid = fopen('anal_3.txt','w');
elseif npt == 'a',
    fid = fopen('anal_3.txt','a');
else
    error('Unexpected answer for input.')
end;

v = zeros(1,sf);

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
    if ~files(o).isdir
        fname = files(o).name;
        ffnm = [where1 fname];
        load(ffnm)
        
% Computing the minimums and sorting the files
        ibiv = IntraBurstIvVar(3);
        if ibiv > n
            eval(['copyfile(''',ffnm,''',''',respath2,''');']);
        else
            eval(['copyfile(''',ffnm,''',''',respath1,''');']);
        end
        
        fprintf(fid,'%s',num2str(ibiv));
        fprintf(fid,'\n');
        
        v(o) = ibiv;
        
    end
    waitbar(o/sf)   %Progress indicator
end

fclose(fid);

close(wb);   %Close progress indicator
cd(mm);