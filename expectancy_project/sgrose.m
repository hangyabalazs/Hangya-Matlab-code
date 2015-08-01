function sgrose
%SGROSE   Rose diagram for delta phase values.
%   SGROSE draws circular histograms ('rose diagram') and phase histograms
%   for delta phase values.
%
%   See also SGPHASERTCORR4.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];     % input directory
resdir = [DATAPATH 'Expectancy\FzCzPz_SE\rose2\'];    % output directory
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'FiveChannel3_phase_rt.mat'];
load(ff)
data = FiveChannel3_phase_rt;    % dimensions: indiv.; cond.; EEG chan.;
                                 % delta phase-sorted delta phase/theta phase/
                                 % reaction time; trials

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
H1 = figure;
H2 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1,H2)
        
%         saveas(H1,['rose_delta ' titlestr '.fig'])    % save
%         saveas(H2,['phasehist_delta ' titlestr '.fig'])
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1,H2)

deltaphase = cell(1,13);
for k = 1:13
    deltaphase{k} = squeeze(data(k,dim2,dim3,1,:))';    % theta phase
end

figure(H1);     % rose diagram
alldelta = cell2mat(deltaphase);
[tt,rr] = rose(ones(1,8));
R = polar(tt,rr/100);
hold on
[tout, rout] = rose(alldelta);
polar(tout,rout/sum(rout))
hold off
delete(R)
title(titlestr)

edges2 = -pi:2*pi/18:pi;    % phase histogram
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
figure(H2);
for k = 1:13
    fr1_indiv = squeeze(data(k,dim2,dim3,1,:));
    fr1_indiv = fr1_indiv(:);
    [nm xout] = histc(fr1_indiv,edges2);
    nm = nm(1:end-1);
    nm = nm / sum(nm);
    plot([cnts2 cnts2+2*pi],[nm' nm'],'Color',[200 200 200]/255,'LineWidth',1)
    hold on
end
fr1_all = squeeze(data(:,dim2,dim3,1,:));
fr1_all = fr1_all(:);
[nm xout] = histc(fr1_all,edges2);
nm = nm(1:end-1);
nm = nm / sum(nm);      % show probabilities
plot([cnts2 cnts2+2*pi],[nm' nm'],'k','LineWidth',2)
ylim([0 0.3])
xlim([min(cnts2) max(cnts2+2*pi)])
hold off
    