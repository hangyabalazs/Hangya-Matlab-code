function b_zshift_xlsconvert
%ZSHIFT_XLSCONVERT   Converts excel output of ZSHIFTRUN.
%   ZSHIFT_XLSCONVERT reads excel output of ZSHIFTRUN and saves a new excel
%   file to a specific location containing the longest segments of cells
%   only.
%
%   See also ZSHIFTRUN2, ZSHIFTRUN3, ZSHIFTRUN_RAO, ZSHIFT_XLSPROCESS
%   and ZSHIFT_XLSCONVERT2.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Zshift_rao\text\'];
fn = [dr 'zshift_theta.xls'];
[mtx ntx atx] = xlsread(fn);

% Select maximal length segment
nmlist = ntx(1:end,1);
cell_list = unique(nmlist);
for i = 1:length(cell_list)
    actinx = find(strcmp(cell_list(i),nmlist));
    len = zeros(1,length(actinx));
    for j = 1:length(actinx)
        len(j) = mtx(actinx(j),2) - mtx(actinx(j),1);
    end
    actinxinx = find(len==max(len));
    inx = actinx(actinxinx);
    atx2{i,1} = ntx{inx,1};
    for k = 2:7
        atx2{i,k} = mtx(inx,k-1);
    end
end

% Write excel file
mm = pwd;
cd(dr)
xlswrite('zshift2_theta_temp',atx2,'maxlen','A2');
cd(mm)