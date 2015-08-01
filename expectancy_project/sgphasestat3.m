function sgphasestat3
%SGPHASESTAT3   Comparison of phase distributions on expectancy data.
%   SGPHASESTAT3 compares delta phases recorded on Fz, Cz and Pz electrodes
%   under different conditions (10%, 37%, 64% and 91%) by Watson-test.
%   Concentration parameters are compared by an F-test due to Fisher (see
%   KAPPACOMPARE for details). Results are plotted and saved along with
%   circular means. Edit code to modify input and output directories!
%
%   See also KAPPACOMPARE and WATSONTWO.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\kappacompare\FIR1_half3\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'delta05-3Hz_theta_singlephase_RT_FIR1_2048.mat'];
load(ff)
data = singleEEGHILB;

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
dim4str = {'delta'};
dim4 = 1;    % delta
H1 = figure;
kappastat = struct('Fz',[],'Cz',[],'Pz',[]);
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim4str{dim4} ' ' dim3str{dim3} ' ' dim2str{dim2}];
        bluestr = dim4str{dim4};
        kappastat = main(data,dim2,dim3,dim4,dim3str,kappastat,titlestr,bluestr,H1);
    end
end

% Save
saveas(H1,'delta_FIR1_half3_kappacomp.fig')
save kappastat2 kappastat
cd(mm)

% -------------------------------------------------------------------------
function kappastat = main(data,dim2,dim3,dim4,dim3str,kappastat,titlestr,bluestr,H1)

% Plot summary
edges2 = -pi:2*pi/18:pi;
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
figure(H1);
subplot(3,4,(dim3-1)*4+dim2)
hold on
fr1_all = squeeze(data(:,dim2,dim3,dim4,:));
fr1_all = fr1_all(:);
[nm xout] = histc(fr1_all,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm' nm'],'b')
ylim([0 200])

% Circular statistics
if dim2 < 4
    fr2_all = squeeze(data(:,dim2+1,dim3,dim4,:));
    fr2_all = fr2_all(:);
    [Fp FH] = kappacompare2(fr2_all,fr1_all,'rad','onesided');
    eval(['kappastat.' dim3str{dim3} '(dim2) = Fp;']);
end
ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment
mn = angle(ftm);   % mean angle
mn = mod(mn,2*pi) * 180 / pi;
mvl = abs(ftm);     % mean resultant length
[kappa_ML kappa] = kappaest(fr1_all,'rad');
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
    [bluestr ' mean: ' num2str(mn)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
    [bluestr ' mvl: ' num2str(mvl)],'Color','blue')
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
    [bluestr ' kappa: ' num2str(kappa)],'Color','blue')
title(titlestr)