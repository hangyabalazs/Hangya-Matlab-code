function aclust_call
%ACLUST_CALL   Caller function for ACLUST.
%   ACLUST_CALL prompts the user to select file, loads discriminated unit,
%   calls cluster analysis and saves results. Edit the program code to
%   modify the input and result directories.
%
%   See also ACLUST and ACLUSTRUN.

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'LentiFenti\disc\'];  % input directory
resdir1 = [DATAPATH 'LentiFenti\Cluster\mat\'];  % result directory
resdir2 = [DATAPATH 'LentiFenti\Cluster\fig\'];
mm = pwd;
cd(inpdir)

% Cluster analysis
[fname pathname] = uigetfile(inpdir);
load(fullfile(pathname,fname))
[Burstiness,IntraBurstFrequency,IntraBurstSpikeNumber,BurstLength,BurstFrequency,...
    Silence,Gap,After,IsiMatrix,Burst,H] = aclust(vdisc);

ach = allchild(H);     % figure title
ax = findobj(ach,'type','axes');
cmps = strread(fname,'%s','delimiter','_');
titlestr = [];
for tt = 1:length(cmps)
    titlestr = [titlestr ' ' cmps{tt}];
end
title(ax(end),titlestr)

% Save
cd(resdir1)
fn = [fname(1:end-6) '_CLUST.mat'];
save(fn,'Burstiness','IntraBurstFrequency','IntraBurstSpikeNumber','BurstLength',...
    'BurstFrequency','Silence','Gap','After','IsiMatrix','Burst','vdisc')
cd(resdir2)
fn2 = [fname(1:end-6) '_CLUST.fig'];
saveas(H,fn2)
close(H)
cd(inpdir)
cd(mm)