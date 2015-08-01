function b_theta(f,newstep,wdir)
%THETA   Theta segment selector.
%   THETA(F,NEWSTEP) loads THETASELECTORRUN output matrices (OM) and saves 'ThetaSegments'
%   matrix, which is a 2-N matrix containing first points of theta segments in its first
%   row and last points of theta segments in its second row. F should be the scalevector 
%   used in WAVELET and NEWSTEP should be the constant which characterizes downsampling 
%   of the EEG (sampled on 10 kHz). 
%
%   You are able to modify the input and the result directory through editing the program
%   code, or you can use THETA(F,NEWSTEP,WDIR) syntax, where WDIR is the working directory
%   of THETA function. WDIR should contain thetaselection data directory - usually 
%   THETASELECTORRUN does this job.
%
%   Theta band: 2.5 - 6 Hz.
%
%   See also THETASELECTORRUN, NONTHETA and THETA_LONG.

% Input arguments check
error(nargchk(2,3,nargin));

% Define directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end
if nargin < 3
    wdir = [DATAPATH,'HCN\Wavelet2\'];
end
where = [wdir '\thetaselection\matrix\'];    %Here are the data files
cd(wdir)
if ~b_isdir2('theta_segments')
    mkdir theta_segments
end
cd('theta_segments');  %Here are the results

% Import
files = dir(where);
files = files(3:end);
sf = length(files);
mm = pwd;

wb = waitbar(0,'Running THETA...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Loading THETASELECTORRUN output matrix
for o = 1:sf
    filenam = files(o).name;
    ff = [where filenam];
    load(ff)
    maxloc = Out(1,:);
    
% Finding theta segments
    fnd = find(f>6);
    pwind1 = fnd(end);
    fnd = find(f<2.5);
    pwind2 = fnd(1);
    fml = find(maxloc>pwind1&maxloc<pwind2);
    if ~isempty(fml)
        th = segm(fml,newstep);
    end
    
% Save
    ThetaSegments = th;
    th = [];
    fs = findstr(filenam,'_');
    fln1 = filenam((fs(2)+1):(fs(3)-1));
    fln2 = filenam((fs(3)+1):end-4);
    eval(['save(''THETA_SEGMENTS_',fln1,'_',fln2,'.mat'',''ThetaSegments'')']);
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator

% ----------------------------------------------------------------------------------------
function th = segm(ip,newstep)
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
preth = zeros(2,lenfdr+1);
preth(1,1) = ip(1);
for t = 1:lenfdr
    preth(2,t) = ip(fdr(t));
    preth(1,t+1) = ip(fdr(t)+1);
end
preth(2,end) = ip(end);
th = preth .* newstep;