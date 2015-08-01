function rim2 = icleanshift(rim)
%ICLEANSHIFT   Deletes sum-vectors from a graph.
%   RIM2 = ICLEANSHIFT(RIM) needs an input argument of a series of
%   adjacency matrices in a 3D array, e.g. time-dependent vectors between
%   fixed points. Non-existent vectors should be represented by zeros or
%   NaNs and real non-zero values should stand for existing vectors.
%   ICLEANSHIFT cleans those vectors from the adjacency matrices that can
%   be generated as the sum of two other vectors. This property is
%   equivalent with the existance of a one value in the raw and a one value
%   in the coloumn with the same (cloloumn or raw) index of the given
%   vector in the adjacency matrix.
%
%   See also ICLEANSHIFT.

% Clean vectors
s3 = size(rim,3);
rim2 = rim;
for k = 1:s3
    am = double(nan2zero(rim(:,:,k))>0);        % adjacency matrix
    sa = size(am);
    nz = find(am);
    for t = 1:length(nz)        % cycle on non-zero elements
        [x,y] = ind2sub(sa,nz(t));
        rw = am(x,:);       % raw of the element
        cl = am(:,y);       % coloumn of the element
        rw(y) = 0;          % replace the element with zero
        cl(x) = 0;
        dp = rw * cl;       % dot product!
        if logical(dp)      % if the dot product is different from zero,
            rim2(x,y,k) = NaN; % a one with the same index in the raw and coloumn exist
        end
    end
end