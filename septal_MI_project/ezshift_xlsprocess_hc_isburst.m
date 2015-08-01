function ezshift_xlsprocess_hc_isburst
%EZSHIFT_XLSPROCESS_HC_ISBURST    Counts bursting neurons.
%   EZSHIFT_XLSPROCESS_HC_ISBURST works on the output excel file of
%   EZSHIFT_XLSCONVERT_HC. It correlates bursting properties with z-shift
%   results.
%
%   See also EZSHIFT_XLSPROCESS.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel files
global DATAPATH
dr = [DATAPATH 'Ezshift_hc\'];
fn = [dr 'zshift_longest_theta.xls'];
clusdir_theta = [DATAPATH 'Burst\Cluster\Theta_hc\'];
mm = pwd;
cd(dr)
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'maxlen');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

dr = [DATAPATH '5Uentropy\'];
fn = [dr 'burstinfo_hc.xls'];
cd(dr)
headerrows = 1;
[mtx2 ntx2 atx2] = xlsread(fn,'burstinfo');
ntx2(1:headerrows,:) = [];
atx2(1:headerrows,:) = [];
names = ntx2(:,1);

% Plot I
zshs_bursting = [];
zshs_nonbursting = [];
H1 = figure;
hold on
for i = 1:size(atx,1)
    nm = atx{i,1};
    sc = strcmp(nm,names);
    t = find(sc);
    if  ~isnan(atx{i,4})
        if atx2{t(1),2}     % bursting
            h3 = plot(atx{i,4},atx{i,5},'.','MarkerSize',30,'Color',[0 1 0]);
            zshs_bursting(end+1) = atx{i,4};
        elseif ~atx2{t(1),2}    % non-bursting
            h3 = plot(atx{i,4},atx{i,5},'.','MarkerSize',30,'Color',[0 0.5 0]);
            zshs_nonbursting(end+1) = atx{i,4};
        end
    end
end
xlim([-10000 10000])
cd(mm)
median(zshs_bursting)
median(zshs_nonbursting)