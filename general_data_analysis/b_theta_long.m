function b_theta_long(f,newstep,wdir)
%THETA_LONG   Version of THETA for long registrations.
%   You have to use THETA_LONG for registrations longer than 200 seconds. THETA_LONG
%   runs on the output of THETASELECTORRUN_LONG. THETA_LONG does exactly the same as 
%   THETA, only the directory names and the saved file names differ. See THETA for
%   details.
%
%   THETA_LONG uses an input and a result directory. You are able to modify them through
%   editing the program code, or you can use THETA_LONG(F,NEWSTEP,WDIR) syntax, where 
%   WDIR is the working directory of THETA_LONG function. WDIR should contain thetaselection
%   data directory - usually THETASELECTORRUN does this job.
%
%   Theta band: 2.5 - 6 Hz.
%
%   See also THETASELECTORRUN, THETASELECTORRUN_LONG, NONTHETA, NONTHETA_LONG and THETA.

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
where = [wdir '\thetaselection_long\matrix\'];    %Here are the data files
cd(wdir)
if ~b_isdir2('theta_segments_long')
    mkdir theta_segments_long
end
cd(['\theta_segments_long\']);  %Here are the results

% Import
files = dir(where);
files = files(3:end);
sf = length(files);
mm = pwd;

wb = waitbar(0,'Running THETA_LONG...','Position',[360 250 275 50]);    %Progress indicator
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
    eval(['save(''THETA_SEGMENTS_LONG',filenam(18),'_',filenam(20:25),'.mat'',''ThetaSegments'')']);
    
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