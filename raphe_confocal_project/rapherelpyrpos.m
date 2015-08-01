function rapherelpyrpos
%RAPHERELPYRPOS   Distance of the pyramidal layer from the image edges.
%   RAPHERELPYRPOS calculates the distance of the pyramidal layer from the
%   image edges in image pixels. It saves the minimal distances over all
%   images in the input directory.
%
%   See also RAPHEPYRPOS and RAPHECONFOCAL_RUN.

% Directories
global DATAPATH
inpdir_gfp = 'X:\Zsolt\zs7\zs7g_rost_denzitas\balazs\905\gfp\';
inpdir_pyrpos = [DATAPATH 'Raphe\zs7\pyrlay\'];
resdir = [DATAPATH 'Raphe\zs7\relpyrpos\'];
mm = pwd;
cd(resdir)
dr = dir(inpdir_gfp);
dr = dr(3:end);
sf = length(dr);

% Load images
prepyr = zeros(1,sf);
postpyr = zeros(1,sf);
for o = 1:sf
    fn_gfp = [inpdir_gfp dr(o).name];
    [Igreen,cmap] = imread(fn_gfp,'tif');
    
% Load pyramidal layer position
    [pth fnm xtn] = fileparts(fn_gfp);
    fn_pp = [inpdir_pyrpos fnm '_PP.mat'];
    load(fn_pp)
    
% Distances
    prepyr(o) = pyrpos - 1;
    postpyr(o) = size(Igreen,2) - pyrpos;
end

% Minimal distances
mprepyr = min(prepyr);
mpostpyr = min(postpyr);

% Save
[pth fnm xtn] = fileparts(fn_gfp);
ff = [fnm(1:8) '_RPP.mat']    % save 'pyrpos'
save(ff,'mprepyr','mpostpyr')
cd(mm)