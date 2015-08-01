function raphepyrlay
%RAPHEPYRLAY   Determine the position of the pyramidal layer.
%   RAPHEPYRLAY loads and displays conocal microscope images. The user is
%   asked to mark the position of tha pyramidal layer by pressing the left
%   mouse button over the center of the layer. The program saves the
%   determined positions.
%
%   See also RAPHECONFOCAL_RUN.

% Directories
global DATAPATH
inpdir_gfp = 'X:\Zsolt\zs7\zs7g_rost_denzitas\balazs\905\gfp\';
resdir = [DATAPATH 'Raphe\zs7\pyrlay\'];
mm = pwd;
cd(resdir)
dr = dir(inpdir_gfp);
dr = dr(3:end);
sf = length(dr);

% Load images
for o = 1:sf
    fn_gfp = [inpdir_gfp dr(o).name];
    [Igreen,cmap] = imread(fn_gfp,'tif');
    Igreen = double(Igreen);
    Igreen = Igreen / max(Igreen(:));
    
% Get pyr. layer position
    H = figure;
    imagesc(Igreen)
    set(gcf,'Position',[45 140 1360 680])
    bp = 1;
    while bp
        disp('Left click over the center of the pyramidal layer!')
        bp = waitforbuttonpress;
        xyco = get(gca,'CurrentPoint');
        hold on
        plot(xyco(1),xyco(3),'w.','MarkerSize',20)
        pause(1)
    end
    
    pyrpos = xyco(1);
    [pth fnm xtn] = fileparts(fn_gfp);
    ff = [fnm '_PP.mat']    % save 'pyrpos'
    save(ff,'pyrpos')
    close(H)
end
cd(mm)