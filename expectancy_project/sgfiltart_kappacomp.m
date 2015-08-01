function sgfiltart_kappacomp
%SGFILTART_KAPPACOMP   Comparison of phase distributions on expectancy data.
%   SGFILTART_KAPPACOMP performs tests for concentration parameter
%   comparison between real and simulated data (see rebuttal to J Neurosci,
%   expectancy project).
%
%   See also SGPHASESTAT3.

% Directories
global DATADIR
global DATAPATH
inpdir1 = [DATADIR 'human_SG\setdataconcatFZCZPZC3C4filtart_ms1400_rand\'];
inpdir2 = [DATADIR 'human_SG\setdataconcatFZCZPZC3C4filtart_ms1400\'];
resdir = [DATAPATH 'Expectancy\filtart_ms1400_rand\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir1 'filtartdelta.mat'];
load(ff)
data1 = erpphase;
ff = [inpdir2 'filtartdelta.mat'];
load(ff)
data2 = erpphase;

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'C3' 'Cz' 'C4' 'Pz'};
kappastat = struct('Fz',[],'C3',[],'Cz',[],'C4',[],'Pz',[]);
for dim2 = 1:4
    for dim3 = 1:5
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        bluestr = 'delta';
        kappastat = main(data1,data2,dim2,dim3,dim3str,kappastat,titlestr,bluestr);
    end
end

% Save
save kappastatnew kappastat
cd(mm)

% -------------------------------------------------------------------------
function kappastat = main(data1,data2,dim2,dim3,dim3str,kappastat,titlestr,bluestr)

% Plot summary
% edges2 = -pi:2*pi/18:pi;
% cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
% figure(H1);
% subplot(length(dim3str),4,(dim3-1)*4+dim2)
% hold on
fr1_all = squeeze(data1(:,dim2,dim3,:));
fr1_all = fr1_all(:);
fr2_all = squeeze(data2(:,dim2,dim3,:));
fr2_all = fr2_all(:);
% nm = histc(fr1_all,edges2);
% nm = nm(1:end-1);
% plot([cnts2 cnts2+2*pi],[nm' nm'],'b')
% ylim([0 200])

% Circular statistics
[Fp FH] = kappacompare2(fr2_all,fr1_all,'rad','onesided');
eval(['kappastat.' dim3str{dim3} '(dim2) = Fp;']);
% ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment
% mn = angle(ftm);   % mean angle
% mn = mod(mn,2*pi) * 180 / pi;
% mvl = abs(ftm);     % mean resultant length
% [kappa_ML kappa] = kappaest(fr1_all,'rad');
% n = length(fr1_all);
% z = n * (mvl ^ 2);  % Rayleigh's Z statistic
% p = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
%     (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
% x_lim = xlim;
% y_lim = ylim;
% text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
%     [bluestr ' mean: ' num2str(mn)],'Color','blue')
% text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,...
%     [bluestr ' mvl: ' num2str(mvl)],'Color','blue')
% text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,...
%     [bluestr ' kappa: ' num2str(kappa)],'Color','blue')
% text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.6,...
%     ['Rayleigh p: ' num2str(p)],'Color','blue')
% title(titlestr)