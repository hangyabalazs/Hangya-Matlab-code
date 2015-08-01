function ezshift_xlsprocess_isburst
%EZSHIFT_XLSPROCESS_ISBURST    Counts bursting neurons.
%   EZSHIFT_XLSPROCESS_ISBURST works on the output excel file of
%   EZSHIFT_XLSCONVERT. It counts bursting cells among MS neurons with
%   Z-shift in/out of the range 0-200 ms.
%
%   See also EZSHIFT_XLSPROCESS.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel files
global DATAPATH
dr = [DATAPATH 'Ezshift\'];
fn = [dr 'zshift_longest_theta.xls'];
mm = pwd;
cd(dr)
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'maxlen');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

dr = [DATAPATH 'Zshift3\text\'];
fn = [dr 'zshift_theta4.xls'];
cd(dr)
headerrows = 1;
[mtx2 ntx2 atx2] = xlsread(fn,'burstinfo');
ntx2(1:headerrows,:) = [];
atx2(1:headerrows,:) = [];
names = ntx2(:,1);

% Plot I
isburst_in = 0;
n_in = 0;
isburst_under = 0;
n_under = 0;
isburst_over = 0;
n_over = 0;
for i = 1:size(atx,1)
    nm = atx{i,1};
    sc = strcmp(nm,names);
    t = find(sc);
    if ~isnan(atx{i,4})
        if atx{i,4} > 0 && atx{i,4} < 2000
            isburst_in = isburst_in + atx2{t(1),2};
            n_in = n_in + 1;
        elseif atx{i,4} <= 0
            isburst_under = isburst_under + atx2{t(1),2};
            n_under = n_under + 1;
        elseif atx{i,4} >= 2000
            isburst_over = isburst_over + atx2{t(1),2};
            n_over = n_over + 1;
        end
    end
end
cd(mm)
isburst_in
isburst_under
isburst_over