function icleanshift2(rim,riml)
%ICLEANSHIFT2   Deletes sum-vectors from a graph.
%   ICLEANSHIFT2(RIM,RIML) performs the same calculation as ICLEANSHIFT on
%   rIMax (RIM) and RIMaxLoc (RIML) matrices. See IMISHIFT for details on
%   input arguments and ICLEANSHIFT for details on calculation. Resulting
%   matrices are saved in a result directory.
%
%   See also ICLEANSHIFT and IMISHIFT.

% Directories
global DATAPATH
resdir = [DATAPATH 'Ulbert\MImap2\'];
mm = pwd;

% Clean vectors
s3 = size(rim,3);
rIMax = rim;
rIMaxLoc = riml;
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
            rIMax(x,y,k) = NaN; % a one with the same index in the raw and coloumn exist
            rIMaxLoc(x,y,k) = NaN;
        end
    end
end

% Save
cd(resdir)
save('MIshiftfine_clean','rIMax','rIMaxLoc')
cd(mm)