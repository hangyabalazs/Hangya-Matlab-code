function sglatrtcorr2
%SGLATRTCORR2   Correlation between reaction time and P300 latency.
%   SGLATRTCORR2 displays and saves  a scatter plot showing reaction time
%   vs. P300 latency.
%
%   See also SGPHASERTCORR4 and SGPHASELATCORR2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\latrt\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'dphase_drt_dlat.mat'];    % dimensions: indiv.; cond.; EEG chan.;
                                        % delta phase/reaction time/P300
                                        % latency; trials
load(ff)
data = dphase_drt_dlat;

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
H1 = figure;
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        main(data,dim2,dim3,titlestr,H1)
        
        saveas(H1,['rt_vs_lat ' titlestr ' FIR1_half3_47.fig'])    % save
    end
end

cd(mm)

% -------------------------------------------------------------------------
function main(data,dim2,dim3,titlestr,H1)

lat = cell(1,13);    % extract rt. and lat.
spvr_delta = zeros(13,100,2);
for k = 1:13
    lat{k} = squeeze(data(k,dim2,dim3,3,:));
    rt_vs_lat = squeeze(data(k,dim2,dim3,[2 3],:))';
    spvr_delta(k,:,1:2) = sortrows(rt_vs_lat,1);
end

figure(H1);     % scatter plot
spvr = spvr_delta;
psphase = spvr(:,:,1);
sphase = psphase(:);
pslat = spvr(:,:,2);
slat = pslat(:);
psp = [sphase slat];
sp = sortrows(psp,1);
sphase = sp(:,1);
slat = sp(:,2);
plot(sphase,slat,'.')
title(titlestr)