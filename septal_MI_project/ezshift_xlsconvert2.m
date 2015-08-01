function ezshift_xlsconvert2
%EZSHIFT_XLSCONVERT2   Converts excel output of EZSHIFTRUN.
%   EZSHIFT_XLSCONVERT2 reads excel output of EZSHIFTRUN and saves a new
%   excel file to a specific location containing the segment with the
%   highest Z-value for each cell.
%
%   See also EZSHIFTRUN, EZSHIFT_XLSPROCESS and EZSHIFT_XLSCONVERT2.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Ezshift\'];
fn = [dr 'zshift_theta.xls'];
headerrows = 1;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Select maximal Z-value segment
nmlist = ntx(1:end,1);
cell_list = unique(nmlist);
for i = 1:length(cell_list)
    actinx = find(strcmp(cell_list(i),nmlist));
    Zmax = zeros(1,length(actinx));
    for j = 1:length(actinx)
        Zmax(j) = mtx(actinx(j),4);
    end
    actinxinx = find(Zmax==max(Zmax));
    inx = actinx(actinxinx);
    atx2{i,1} = ntx{inx,1};
    for k = 2:9
        if k < 8
            atx2{i,k} = mtx(inx,k-1);
        else
            atx2{i,k} = atx{inx,k};
        end
    end
end

% Write excel file
mm = pwd;
cd(dr)
xlswrite('zshift_maxz_theta',atx2,'maxz','A2');
cd(mm)