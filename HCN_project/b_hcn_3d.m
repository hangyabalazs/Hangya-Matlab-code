function b_hcn_3d
%HCN_3D   3D plot from burst analysis results.
%   HCN_3D plots one point per analysed segment in a 3D coordinate system.
%   X axis: minimum of first-spike-CV values from 2 to 7 clusters
%   Y axis: burst length CV
%   Z axis: intraburst interval CV
%
%   See also ITCLUST2.

% Input argument check
error(nargchk(0,0,nargin))

% Input directory
cd d:\_analysis\matlab_data\HCN\Analysis2

% 3D plot
d = dir;
for i = 3:length(d)
    if ~d(i).isdir
        if ~isempty(findstr(d(i).name,'ICA'))
            load(d(i).name)
            
            mfsp = min(FirstSpikeCv(2:7));
            im = find(FirstSpikeCv==mfsp);
            bl = BurstLengthCv(im(1));
            ib = IntraBurstIvCv(im(1));
            
            plot3(mfsp,bl,ib,'.','MarkerSize',10);
            hold on
            grid on
        end
    end
end