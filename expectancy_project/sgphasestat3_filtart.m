function sgphasestat3_filtart
%SGPHASESTAT3_FILTART   Comparison of phase distributions on expectancy data.
%   SGPHASESTAT3_FILTART is a version of SGPHASESTAT3 on data created
%   to test possible filter artifacts in Experiment1. All details can be
%   found in the documentation of SGPHASESTAT3.
%
%   See also SGPHASESTAT3.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\setdataconcatFZCZPZC3C4nofilt_halfepoch\'];
resdir = [DATAPATH 'Expectancy\filtart_nofilt_halfepoch\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'filtartdelta.mat'];
load(ff)
% data = singleEEGHILB;
data = erpphase;

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'C3' 'Cz' 'C4' 'Pz'};
kappastat = struct('Fz',[],'C3',[],'Cz',[],'C4',[],'Pz',[]);
H1 = figure;
for dim2 = 1:4
    for dim3 = 1:5
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        bluestr = 'delta';
        kappastat = main(data,dim2,dim3,dim3str,kappastat,titlestr,bluestr,H1);
    end
end

% Save
saveas(H1,'delta_kappacomp.fig')
save kappastat2 kappastat
cd(mm)

% -------------------------------------------------------------------------
function kappastat = main(data,dim2,dim3,dim3str,kappastat,titlestr,bluestr,H1)

% Plot summary
edges2 = -pi:2*pi/18:pi;
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
figure(H1);
subplot(length(dim3str),4,(dim3-1)*4+dim2)
hold on
fr1_all = squeeze(data(:,dim2,dim3,:));
fr1_all = fr1_all(:);
nm = histc(fr1_all,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm' nm'],'b')
ylim([0 200])

% Circular statistics
if dim2 < 4
    fr2_all = squeeze(data(:,dim2+1,dim3,:));
    fr2_all = fr2_all(:);
    [Fp FH] = kappacompare2(fr2_all,fr1_all,'rad','onesided');
    eval(['kappastat.' dim3str{dim3} '(dim2) = Fp;']);
end
ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment
mn = angle(ftm);   % mean angle
mn = mod(mn,2*pi) * 180 / pi;
mvl = abs(ftm);     % mean resultant length
[kappa_ML kappa] = kappaest(fr1_all,'rad');
n = length(fr1_all);
z = n * (mvl ^ 2);  % Rayleigh's Z statistic
p = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
    (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    [bluestr ' mean: ' num2str(mn)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    [bluestr ' mvl: ' num2str(mvl)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    [bluestr ' kappa: ' num2str(kappa)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.6,...
    ['Rayleigh p: ' num2str(p)],'Color','blue')
title(titlestr)