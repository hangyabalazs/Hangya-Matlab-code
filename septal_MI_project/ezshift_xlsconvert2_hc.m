function ezshift_xlsconvert2_hc
%EZSHIFT_XLSCONVERT2_HC   Converts excel output of EZSHIFTRUN_HC.
%   EZSHIFT_XLSCONVERT2_HC reads excel output of EZSHIFTRUN_HC and saves a
%   new excel file to a specific location containing the segment with the
%   highest Z-value for each cell.
%
%   See also EZSHIFTRUN_HC, EZSHIFT_XLSPROCESS_HC and
%   EZSHIFT_XLSCONVERT2_HC.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Ezshift_hc\'];
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
    for k = 2:7
        atx2{i,k} = mtx(inx,k-1);
    end
end

% Write excel file
mm = pwd;
cd(dr)
xlswrite('zshift_maxz_theta',atx2,'maxz','A2');
cd(mm)