function sglatrtcorr
%SGLATRTCORR   Correlation between reaction time and P300 latency.
%   SGLATRTCORR displays and saves  a scatter plot showing reaction time
%   vs. P300 latency.
%
%   See also SGPHASERTCORR4 and SGPHASELATCORR2.

% Stucture of singleEEGHILB:
% size:
%    13     4     3     3   100
% dim1 = 13 subjects
% dim2 = 4 conditions, 1 -10%, 2 -37%, 3 -64% 4 -91%
% dim3 = 3 EEG channels, 1-Fz, 2-Cz, 3-Pz
% dim4= 1: delta phase-sorted delta phase, 2: delta phase-sorted theta 
%   phase, 3: delta phase-sorted RT
% dim5 = 100 trials

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\latrt\FIR1_half3_47\'];
mm = pwd;
cd(resdir)

% Import
ff = [inpdir 'phase_amp_lat_RT_FIR1_2048_12_08_2008.mat'];
load(ff)
data = singP3lathilbsort;
data2 = avgRT;

% Main
dim2str = {'10%' '37%' '64%' '91%'};
dim3str = {'Fz' 'Cz' 'Pz'};
plat = zeros(13,100);
prt = zeros(13,100);
for dim2 = 1:4
    for dim3 = 1:3
        titlestr = [dim3str{dim3} ' ' dim2str{dim2}];
        for k = 1:13
            plat(k,:) = squeeze(data(k,dim2,dim3,:))';
            prt(k,:) = squeeze(data2(k,dim2,:))';
        end
        rt = prt(:);
        lat = plat(:);
        H1 = plot(rt,lat,'.');
        saveas(H1,['rt_vs_lat ' titlestr ' FIR1_half3_47.fig'])    % save
    end
end

cd(mm)