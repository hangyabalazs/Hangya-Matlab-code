function ezshift_xlsconvert
%EZSHIFT_XLSCONVERT   Converts excel output of EZSHIFTRUN.
%   EZSHIFT_XLSCONVERT reads excel output of EZSHIFTRUN and saves a new
%   excel file to a specific location containing the longest segments of
%   cells only.
%
%   See also EZSHIFTRUN, EZSHIFT_XLSPROCESS and EZSHIFT_XLSCONVERT2.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Ezshift_Rao\'];
fn = [dr 'zshift_theta.xls'];
headerrows = 1;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

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
xlswrite('zshift_longest_theta',atx2,'maxlen','A2');
cd(mm)