function sgphasestat2
%SGPHASESTAT2   Comparison of phase distributions on expectancy data.
%   SGPHASESTAT2 compares theta, delta and alpha phases recorded on Fz, Cz
%   and Pz electrodes under different conditions (10%, 37%, 64% and 91%) by
%   Watson-test. Results are plotted and saved along with circular means.
%   Edit code to modify input and output directories!
%
%   See also CIRCULAR_MEAN3 and WATSONTWO.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\phasediff\bpF_13_47\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'delta1-3Hz_theta4-7Hz_singlephase_badpassFilter.mat'];
load(ff)
data = singleEEGHILB;

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
dim4str = {'delta' 'theta'};
dim4 = [2 1];    % theta-delta
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim4str{dim4(1)} '-' dim4str{dim4(2)} ' ' dim3str{dim3} ' ' ...
            dim2str{dim2}]
        redstr = dim4str{dim4(1)};
        bluestr = dim4str{dim4(2)};
        main(data,dim2,dim3,dim4,titlestr,redstr,bluestr)
    end
end
cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,dim4,titlestr,redstr,bluestr)

% Plot individuals
fr1 = cell(1,13);
edges = -pi:2*pi/9:pi;
cnts = (edges(1:end-1) + edges(2:end)) / 2;
H1 = figure;
hold on
for k = 1:13
    fr1{k} = squeeze(data(k,dim2,dim3,dim4(1),:));
    [nm xout] = histc(fr1{k},edges);
    nm = nm(1:end-1);
    plot([cnts cnts+2*pi],[nm' nm'],'r')
end
fr2 = cell(1,13);
for k = 1:13
    fr2{k} = squeeze(data(k,dim2,dim3,dim4(2),:));
    [nm xout] = histc(fr2{k},edges);
    nm = nm(1:end-1);
    plot([cnts cnts+2*pi],[nm' nm'],'b')
end
title(titlestr)

% Plot summary
edges2 = -pi:2*pi/18:pi;
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
H2 = figure;
hold on
fr1_all = squeeze(data(:,dim2,dim3,dim4(1),:));
fr1_all = fr1_all(:);
[nm xout] = histc(fr1_all,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm' nm'],'r')
fr2_all = squeeze(data(:,dim2,dim3,dim4(2),:));
fr2_all = fr2_all(:);
[nm xout] = histc(fr2_all,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm' nm'],'b')

% Circular statistics
[U2 p] = b_watsontwo(fr1_all,fr2_all);
mn1 = b_circular_mean3(fr1_all) * 180 / pi;
mn2 = b_circular_mean3(fr2_all) * 180 / pi;
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    [redstr ' mean: ' num2str(mn1)],'Color','red')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    [bluestr ' mean: ' num2str(mn2)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    ['p < ' num2str(p(2))],'Color','black')
title(titlestr)

% Save
saveas(H1,[titlestr '.fig'])
saveas(H2,[titlestr ' summary.fig'])