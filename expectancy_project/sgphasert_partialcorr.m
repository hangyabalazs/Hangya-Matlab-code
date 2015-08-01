function sgphasert_partialcorr
%SGPHASERT_PARTIALCORR   Partial correlation between delta phase and reaction time.
%   SGPHASERT_PARTIALCORR calculates linear-circular partial correlation 
%   between 0 ms delta phase and reaction time (rt.) conditioned on P300 
%   amplitude/latency. Correlation coefficients and p values for
%   significance testing are saved in an Excel table.
%
%   See also SGPHASERTCORR4.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];     % input directory
resdir = [DATAPATH 'Expectancy\phase_rt_partial\'];    % output directory
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'FiveChannel3_phase_rt.mat'];
load(ff)
data = FiveChannel3_phase_rt;    % dimensions: indiv.; cond.; EEG chan.;
                                 % delta phase-sorted delta phase/theta phase/
                                 % reaction time; trials
ff = [inpdir 'FiveChannel3_phase_amp_lat.mat'];
load(ff)
data_amp = FiveChannel3_phase_amp_lat;  % dimensions: indiv.; cond.; EEG chan.;
                                        % delta phase/P300 amp./P300 latency;
                                        % trials
data2 = zeros(13,4,3,4,100);
data2(:,:,:,1:3,:) = data;
amplat = 'amp';
aorl = 2;   % 2: amp; 3: lat.
data2(:,:,:,4,:) = data_amp(:,:,:,aorl,:);    % dimensions: indiv.; cond.; EEG chan.;
                                 % delta phase-sorted delta phase/theta phase/
                                 % reaction time/P300 amp. or lat.; trials

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
xlsout = [{'channel'} {'condition'} {'partial R'} {'p'} ...
    {'R_rt,phase (circular)'} {'p'} {'R_rt,phase (linear)'} {'p'} ...
    {'R_amp,phase (circular)'} {'p'} {'R_amp,phase (linear)'} {'p'}];
for dim3 = 1:3
    for dim2 = 1:4
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        [rho pval R_rt_phase_circular p_rt_phase_circular ...
            R_rt_phase_linear p_rt_phase_linear ...
            R_amp_phase_circular p_amp_phase_circular ...
            R_amp_phase_linear p_amp_phase_linear] = ...
            main(data2,dim2,dim3,titlestr);
        
%         saveas(H1,['sorted_deltaphase2 ' titlestr '.fig'])    % save
%         saveas(H2,['sorted_thetaphase2 ' titlestr '.fig'])
%         saveas(H3,['sorted_phasediff2 ' titlestr '.fig'])
        xlsout(end+1,:) = [{dim3str{dim3}} {dim2str{dim2}} {rho} {pval} ...
            {R_rt_phase_circular} {p_rt_phase_circular} ...
            {R_rt_phase_linear} {p_rt_phase_linear} ...
            {R_amp_phase_circular} {p_amp_phase_circular} ...
            {R_amp_phase_linear} {p_amp_phase_linear}];
    end
end
xlswrite('partial_corr_phase_rt_' amplat '2',xlsout)
cd(mm)

% -------------------------------------------------------------------------
function [rho pval R_rt_phase_circular p_rt_phase_circular ...
    R_rt_phase_linear p_rt_phase_linear ...
    R_amp_phase_circular p_amp_phase_circular ...
    R_amp_phase_linear p_amp_phase_linear] = main(data,dim2,dim3,titlestr)

deltaphase = cell(1,13);
rt = cell(1,13);
amp = cell(1,13);
spvr_delta = zeros(13,100,2);
spvr_theta = zeros(13,100,2);
spvr_diff = zeros(13,100,2);
for k = 1:13
    deltaphase{k} = squeeze(data(k,dim2,dim3,1,:));    % delta phase
    rt{k} = squeeze(data(k,dim2,dim3,3,:));            % reaction time
    amp{k} = squeeze(data(k,dim2,dim3,4,:));           % P300 amplitude
            
    deltaphase_vs_rt_vs_amp = squeeze(data(k,dim2,dim3,[1 3 4],:))';
    spvr_delta(k,:,1:3) = sortrows(deltaphase_vs_rt_vs_amp,1);
end

spvr = spvr_delta;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psamp = spvr(:,:,3);
samp = psamp(:);
psp = [sphase srt samp];
sp = sortrows(psp,1);
sphase = sp(:,1);
srt = sp(:,2);
samp = sp(:,3);

% figure
[b,bint,r,rint,stats] = regress(srt,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);      % correlation
R_rt_phase_circular = sqrt(stats(1));
p_rt_phase_circular = stats(3);
[b,bint,r,rint,stats] = regress(srt,[ones(length(sphase),1) ...
    sphase]);      % correlation
R_rt_phase_linear = sqrt(stats(1));
p_rt_phase_linear = stats(3);
% subplot(1,2,1)
% plot(sphase,1:length(sphase),'k')
% x_lim = xlim;
% ylim([0 length(sphase)])
% y_lim = ylim;
% pl = oom(p);
% text(x_lim(1)+(x_lim(2)-x_lim(1))*0.1,y_lim(1)+(y_lim(2)-y_lim(1))*0.95,...
%     ['R = ' num2str(R) '  p ' pl],'Color','red')
% subplot(1,2,2)
% [sslat sds] = smooth(srt,'circular');
% lssl = length(sslat);
% plot(sslat,1:lssl,'k')
% hold on
% plot(sslat+sds,1:lssl,'Color',[204 204 204]/255)
% plot(sslat-sds,1:lssl,'Color',[204 204 204]/255)
% hold off
% ylim([0 lssl])
% title(titlestr)
% 
% figure
[b,bint,r,rint,stats] = regress(samp,[ones(length(sphase),1) ...
    sin(sphase) cos(sphase)]);      % correlation
R_amp_phase_circular = sqrt(stats(1));
p_amp_phase_circular = stats(3);
[b,bint,r,rint,stats] = regress(samp,[ones(length(sphase),1) ...
    sphase]);      % correlation
R_amp_phase_linear = sqrt(stats(1));
p_amp_phase_linear = stats(3);
% subplot(1,2,1)
% plot(sphase,1:length(sphase),'k')
% x_lim = xlim;
% ylim([0 length(sphase)])
% y_lim = ylim;
% pl = oom(p);
% text(x_lim(1)+(x_lim(2)-x_lim(1))*0.1,y_lim(1)+(y_lim(2)-y_lim(1))*0.95,...
%     ['R = ' num2str(R) '  p ' pl],'Color','red')
% subplot(1,2,2)
% [sslat sds] = smooth(samp,'circular');
% lssl = length(sslat);
% plot(sslat,1:lssl,'k')
% hold on
% plot(sslat+sds,1:lssl,'Color',[204 204 204]/255)
% plot(sslat-sds,1:lssl,'Color',[204 204 204]/255)
% hold off
% ylim([0 lssl])
% title(titlestr)

[rho pval] = lincirc_parcor(srt,sphase,samp);


% -------------------------------------------------------------------------
function [X2 S] = smooth(X,str)

n = 101;
nn = (n - 1) / 2;
m = length(X);
X2 = zeros(size(X));
S = zeros(size(X));
switch str
    case 'linear'
        for k = 1:m
            ind1 = max(1,k-nn);
            ind2 = min(k+nn,m);
            X2(k) = mean(X(ind1:ind2));
            S(k) = std(X(ind1:ind2));
        end
    case 'circular'
        for k = 1:m
            if k - nn < 1
                X2(k) = mean([X(mod2(k-nn,m):m); X(1:k+nn)]);
                S(k) = std([X(mod2(k-nn,m):m); X(1:k+nn)]) / sqrt(n);
            elseif k + nn > m
                X2(k) = mean([X(k-nn:m); X(1:mod2(k+nn,m))]);
                S(k) = std([X(k-nn:m); X(1:mod2(k+nn,m))]) / sqrt(n);
            else
                X2(k) = mean(X(k-nn:k+nn));
                S(k) = std(X(k-nn:k+nn)) / sqrt(n);
            end
        end
end

% -------------------------------------------------------------------------
function r = mod2(a,m)

r = mod(a,m);
if r == 0
    r = m;
end

% -------------------------------------------------------------------------
function dt = pmpi(dt)

dt(dt<-pi) = dt(dt<-pi) + 2 * pi;
dt(dt>pi) = dt(dt>pi) - 2 * pi;

% -------------------------------------------------------------------------
function pl = oom(p)

if isequal(p,0)
    pl = '< 10^{-15}';
elseif p < 0.05 && p > 0.01
    pl = '< 0.05';
elseif p < 0.01 && p > 0.001
    pl = '< 0.01';
elseif p < 0.001 && p > 0.0001
    pl = '< 0.001';
elseif p < 0.0001
    ep = (-3:-1:-16);
    ords = 10 .^ ep;
    iop = find(ords>p,1,'last');
    epp = ep(iop);
    pl = ['< 10^{' num2str(epp) '}'];
else
    pl = ['= ' num2str(p)];
end 