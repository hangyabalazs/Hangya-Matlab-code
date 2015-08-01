function aclustrun
%ACLUSTRUN    Runs ACLUST on a sequence of files.
%   ACLUSTRUN loads discriminated unit, calls cluster analysis and saves
%   results. Edit the program code to modify the input and result
%   directories.
%
%   See also ACLUST.

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'LentiFenti\disc\'];  % input directory
resdir1 = [DATAPATH 'LentiFenti\Cluster\mat\'];  % result directory
resdir2 = [DATAPATH 'LentiFenti\Cluster\fig\'];
mm = pwd;
cd(inpdir)

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Cluster analysis
for o = 1:sf
    fname = files_short{o};
    load(fname)
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
end
cd(mm)

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
vrs = version;
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);