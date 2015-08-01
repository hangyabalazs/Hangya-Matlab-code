function sgphasertcorr3
%SGPHASERTCORR3   Correlation between delta/theta phase and reaction time.
%   SGPHASERTCORR displayes and saves 0 ms delta, theta and delta minus
%   theta phase histograms with mean reaction time in each phase bin.
%
%   See also SGPHASERTCORR2, SGPHASERTCORR4 and SGPHASERTCORR5.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];     % inpur directory
resdir = [DATAPATH 'Expectancy\C3C42\phase_rt\'];    % output directory
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'FiveChannels_phase_rt.mat'];
load(ff)
data = FiveChannels_phase_rt;    % dimensions: indiv.; cond.; EEG chan.;
                                 % delta phase-sorted delta phase/theta phase/
                                 % reaction time; trials

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'C3' 'Cz' 'C4'};
H1 = figure;
H2 = figure;
H3 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1,H2,H3)
        
        saveas(H1,['hist_deltaphase2 ' titlestr '.fig'])    % save
%         saveas(H2,['hist_thetaphase2 ' titlestr '.fig'])
%         saveas(H3,['hist_phasediff2 ' titlestr '.fig'])
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1,H2,H3)

thetaphase = cell(1,13);
deltaphase = cell(1,13);
phasediff = cell(1,13);
rt = cell(1,13);
spvr_delta = zeros(13,100,2);
spvr_theta = zeros(13,100,2);
spvr_diff = zeros(13,100,2);
for k = 1:13
    deltaphase{k} = squeeze(data(k,dim2,dim3,1,:));     % delta phase
    thetaphase{k} = squeeze(data(k,dim2,dim3,2,:));     % theta phase
    phasediff{k} = deltaphase{k} - thetaphase{k};       % phase difference
    rt{k} = squeeze(data(k,dim2,dim3,3,:));             % reaction time
            
    deltaphase_vs_rt = squeeze(data(k,dim2,dim3,[1 3],:))';
    spvr_delta(k,:,1:2) = sortrows(deltaphase_vs_rt,1);
    thetaphase_vs_rt = squeeze(data(k,dim2,dim3,[2 3],:))';
    spvr_theta(k,:,1:2) = sortrows(thetaphase_vs_rt,1);
    phasediff_vs_rt = [pmpi(deltaphase_vs_rt(:,1)-thetaphase_vs_rt(:,1))...
        thetaphase_vs_rt(:,2)];
    spvr_diff(k,:,1:2) = sortrows(phasediff_vs_rt,1);
end
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers

figure(H1);     % phase hist.
spvr = spvr_delta;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psp = [sphase srt];
spvr = sortrows(psp,1);
[mn_phase mn_rt] = dhist(spvr,edges);
subplot(2,1,1)
plot([cnts cnts+2*pi],[mn_phase mn_phase],'r')
subplot(2,1,2)      % conditional rt.: E(rt|phi1<phase<phi2)
plot([cnts cnts+2*pi],[mn_rt mn_rt],'r')
title(titlestr)

figure(H2);     % phase hist.
spvr = spvr_theta;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psp = [sphase srt];
spvr = sortrows(psp,1);
[mn_phase mn_rt] = dhist(spvr,edges);
subplot(2,1,1)
plot([cnts cnts+2*pi],[mn_phase mn_phase],'r')
subplot(2,1,2)      % conditional rt.: E(rt|phi1<phase<phi2)
plot([cnts cnts+2*pi],[mn_rt mn_rt],'r')
title(titlestr)

figure(H3);     % phase hist.
spvr = spvr_diff;
psphase = spvr(:,:,1);
sphase = psphase(:);
psrt = spvr(:,:,2);
srt = psrt(:);
psp = [sphase srt];
spvr = sortrows(psp,1);
[mn_phase mn_rt] = dhist(spvr,edges);
subplot(2,1,1)
plot([cnts cnts+2*pi],[mn_phase mn_phase],'r')
subplot(2,1,2)      % conditional rt.: E(rt|phi1<phase<phi2)
plot([cnts cnts+2*pi],[mn_rt mn_rt],'r')
title(titlestr)

% -------------------------------------------------------------------------
function dt = pmpi(dt)

dt(dt<-pi) = dt(dt<-pi) + 2 * pi;
dt(dt>pi) = dt(dt>pi) - 2 * pi;

% -------------------------------------------------------------------------
function [mn_phase mn_rt] = dhist(spvr,edges)

n = length(edges);
mn_phase = zeros(1,n-1);
mn_rt = zeros(1,n-1);
for k = 2:n
    inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<edges(k));
    mn_phase(k-1) = length(inx);
    mn_rt(k-1) = mean(spvr(inx,2));
end