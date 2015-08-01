function ezshift_xlsprocess3
%EZSHIFT_XLSPROCESS3    Generates plots from ZSHIFT output excel file.
%   EZSHIFT_XLSPROCESS3 works on the output excel file of EZSHIFTRUN. It
%   generates various plots.
%
%   See also EZSHIFTRUN and EZSHIFT_XLSPROCESS2.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
dr = [DATAPATH 'Ezshift\'];
fn = [dr 'zshift_theta.xls'];
mm = pwd;
cd(dr)
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'Sheet1');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Plot I
PVp = [];
PVn = [];
na = [];
H1 = figure;
hold on
for i = 1:size(atx,1)
    if isequal(atx{i,8},'PVp') & ~isnan(atx{i,4})
        h1 = plot(atx{i,4},atx{i,5},'r.','MarkerSize',50);
        PVp(end+1) = atx{i,4};
    elseif isequal(atx{i,8},'PVn') & ~isnan(atx{i,4})
        h2 = plot(atx{i,4},atx{i,5},'.','MarkerSize',50,'Color',[1 0.69 0.39]);
        PVn(end+1) = atx{i,4};
    elseif isnan(atx{i,8}) & ~isnan(atx{i,4})
        h3 = plot(atx{i,4},atx{i,5},'k.','MarkerSize',30);
        na(end+1) = atx{i,4};
    end
end
saveas(H1,'zshift_allPV.fig')

H2 = mibox(PVp,PVn,'PV+','PV-');
saveas(H2,'zshift_allPVbox.fig')

% Plot II
HCNp = [];
HCNn = [];
na = [];
H3 = figure;
hold on
for i = 1:size(atx,1)
    if isequal(atx{i,9},'HCNp') & ~isnan(atx{i,4})
        h1 = plot(atx{i,4},atx{i,5},'b.','MarkerSize',50);
        HCNp(end+1) = atx{i,4};
    elseif isequal(atx{i,9},'HCNn') & ~isnan(atx{i,4})
        h2 = plot(atx{i,4},atx{i,5},'.','MarkerSize',50,'Color',[51 204 255]/256);
        HCNn(end+1) = atx{i,4};
    elseif isnan(atx{i,9}) & ~isnan(atx{i,4})
        h3 = plot(atx{i,4},atx{i,5},'k.','MarkerSize',30);
        na(end+1) = atx{i,4};
    end
end
saveas(H3,'zshift_allHCN.fig')

H4 = mibox(HCNp,HCNn,'HCN+','HCN-');
saveas(H4,'zshift_allHCNbox.fig')
cd(mm)



% -------------------------------------------------------------------------
function H = mibox(m1,m2,l1,l2)

H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{l1} {l2}]);
[Wp,Wh] = b_ranksum2(m1,m2,'alpha',0.05);
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