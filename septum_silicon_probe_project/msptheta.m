function msptheta(f,newstep,wdir)
%MSPTHETA   Theta segment selector.
%   MSPTHETA(F,NEWSTEP) is modified from THETA (see the help of the latter
%   for details). It calculates non-theta segments as the complementary of
%   theta segments.
%
%   Theta band: 2.5 - 6 Hz.
%
%   See also THETA and MSPTHETA_MAIN2.

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
if ~b_isdir2('nontheta_segments')
    mkdir nontheta_segments
end

% Import
files = dir(where);
files = files(3:end);
sf = length(files);
mm = pwd;

wb = waitbar(0,'Running MSPTHETA...','Position',[360 250 275 50]);    %Progress indicator
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
    cd('theta_segments');
    ThetaSegments = th;
    th = [];
    fs = findstr(filenam,'_');
    fln1 = filenam((fs(2)+1):(fs(3)-1));
    fln2 = filenam((fs(3)+1):end-4);
    eval(['save(''THETA_SEGMENTS_',fln1,'_',fln2,'.mat'',''ThetaSegments'')']);
    cd ..
    
% Finding non-theta segments
    fml = find(~(maxloc>pwind1&maxloc<pwind2));
    if ~isempty(fml)
        nth = segm(fml,newstep);
    end
    
% Save
    cd('nontheta_segments');
    NonThetaSegments = nth;
    nth = [];
    eval(['save(''NONTHETA_SEGMENTS_',fln1,'_',fln2,'.mat'',''NonThetaSegments'')']);
    cd ..
    
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