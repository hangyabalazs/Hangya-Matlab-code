function [Ang, Mvl, Fname] = axls
%AXLS   Collects result in one matrix.
%
%   See also APHASERUN2.

dr = 'D:\_analysis\matlab_data\Hajni\Phase2\';
files = dir(dr);
files = files(3:end);
Ang = [];
Mvl = [];
Fname = {};
for k = 1:length(files)
    if ~files(k).isdir
        fname = files(k).name;
        if isequal(fname(end-3:end),'.mat')
            load([dr fname])
            Ang(end+1) = ang;
            Mvl(end+1) = mvl;
            Fname{end+1} = fname(1:end-10);
        end
    end
end
Ang = Ang';
Mvl = Mvl';
Fname = Fname';