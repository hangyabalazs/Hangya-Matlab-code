function ezshift_xlsprocess_hc
%EZSHIFT_XLSPROCESS_HC    Generates plots from ZSHIFT output excel file.
%   EZSHIFT_XLSPROCESS_HC works on the output excel file of
%   EZSHIFT_XLSCONVERT_HC. It generates various plots.
%
%   See also EZSHIFTRUN_HC and EZSHIFT_XLSPROCESS2_HC.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Ezshift_hc\'];
fn = [dr 'zshift_longest_theta.xls'];
mm = pwd;
cd(dr)
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'maxlen');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Plot I
zshs = [];
zshs2 = [];
H1 = figure;
hold on
for i = 1:size(atx,1)
    h3 = plot(atx{i,4},atx{i,5},'.','MarkerSize',50,'Color',[0 0.5 0]);
    zshs(end+1) = atx{i,4};
    if ~(atx{i,4} < -150 || atx{i,4} > 700)
        zshs2(end+1) = atx{i,4};
    end
end
xlim([-10000 10000])
% saveas(H1,'zshift_maxlen.fig')
% 
% H2 = mibox(zshs,'Z-shift (maxlen)');
% saveas(H2,'zshift_maxlenbox.fig')
cd(mm)



% -------------------------------------------------------------------------
function H = mibox(m1,l1)

H = figure;
boxplot(m1,zeros(size(m1)),'labels',[{l1}]);
[Wp,Wh] = b_ranksum2(m1,zeros(size(m1)),'alpha',0.05);
if Wh
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')