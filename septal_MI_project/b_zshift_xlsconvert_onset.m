function b_zshift_xlsconvert_onset
%ZSHIFT_XLSCONVERT_ONSET   Converts excel output of ZSHIFTRUN_ONSET.
%   ZSHIFT_XLSCONVERT_ONSET reads excel output of ZSHIFTRUN and saves a new
%   excel file to a specific location containing the longest segments of
%   cells only.
%
%   See also ZSHIFTRUN2, ZSHIFTRUN3, ZSHIFTRUN_ONSET, ZSHIFT_XLSPROCESS
%   and ZSHIFT_XLSCONVERT2.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Zshift_onset\text\'];
fn = [dr 'zshift_theta.xls'];
[mtx ntx atx] = xlsread(fn);
dr2 = [DATAPATH 'Zshift4\text\'];
fn2 = [dr2 'zshift2_theta.xls'];
[mtx2 ntx2 atxx2] = xlsread(fn2);
mtxx2 = mtx2(:,1);

% Select maximal length segment
nmlist = ntx(1:end,1);
cell_list = unique(nmlist);
for i = 1:length(cell_list)
    actinx = find(strcmp(cell_list(i),nmlist));
    lena = length(actinx);
    start = zeros(1,lena);
    slg = zeros(1,lena);
    for j = 1:lena
        start(j) = mtx(actinx(j),1);
        lenm = length(mtxx2);
        cmp = zeros(1,lenm);
        for k = 1:lenm
            cmp(k) = isequal(start(j),mtxx2(k));
        end
        slg(j) = ~isempty(find(cmp));
        if slg(j)
            actinxinx = j;
        end
    end
    inx = actinx(actinxinx);
    atx2{i,1} = ntx{inx,1};
    for k = 2:7
        atx2{i,k} = mtx(inx,k-1);
    end
end

% Write excel file
mm = pwd;
cd(dr)
xlswrite('zshift2_theta',atx2,'maxlen','A2');
cd(mm)