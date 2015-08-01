function b_zshift_xlsconvert2
%ZSHIFT_XLSCONVERT   Converts excel output of ZSHIFTRUN.
%   ZSHIFT_XLSCONVERT reads excel output of ZSHIFTRUN and saves a new excel
%   file to a specific location containing the segment with the highest 
%   Z-value for each cell.
%
%   See also ZSHIFTRUN2, ZSHIFTRUN3, ZSHIFTRUN_RAO, ZSHIFT_XLSPROCESS
%   and ZSHIFT_XLSCONVERT.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
fn = [DATAPATH 'Zshift4\text\zshift_theta.xls'];
[mtx ntx atx] = xlsread(fn);

% Select maximal Z-value segment
nmlist = ntx(2:end,1);
cell_list = unique(nmlist);
for i = 1:length(cell_list)
    actinx = find(strcmp(cell_list(i),nmlist));
    Zmax = zeros(1,length(actinx));
    for j = 1:length(actinx)
        Zmax(j) = mtx(actinx(j),4);
    end
    actinxinx = find(Zmax==max(Zmax));
    inx = actinx(actinxinx);
    atx2{i,1} = ntx{inx(1)+1,1};
    for k = 2:9
        atx2{i,k} = mtx(inx(1),k-1);
    end
end

% Write excel file
cd('D:\_analysis\matlab_data\Zshift4\text\')
xlswrite('zshift2_theta',atx2,'maxZ','A2');