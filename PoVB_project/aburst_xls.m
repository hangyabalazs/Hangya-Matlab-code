function aburst_xls
%ABURST_XLS   Converts burst mat files to excel sheets.
%   ABURST_XLS loads burst parameters and saves them to excel files. An
%   xls file is created for each animal.
%
%   See also ACLUSTERCUT.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\control\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\Burst_xls\control\'];
mm = pwd;

% Filelist
[files files_short] = filelist2(inpdir);
sf = length(files_short);

% Main
for o = 1:sf
    fname = files_short{o};     % load
    ff2 = [inpdir fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    cd(resdir)
    cmps = strread(fname,'%s','delimiter','_');
    xlsname = cmps{1};
    if exist([xlsname '.xls'],'file')
        ntx = xlsread(xlsname,'sheet1');
        nty = xlsread(xlsname,'sheet2');
        [ntz mtz atz] = xlsread(xlsname,'sheet3');
    else
        ntx = [];
        nty = [];
        ntz = [];
        atz = {};
        hrow = {[] 'Bness' 'IBF mean' 'IBF sd' 'IBSN mean' 'IBSN sd' 'BL mean' 'BL sd' 'BF'};
        hrowy = {[] 'IBF' 'IBSN' 'BL'};
        xlswrite(xlsname,hrow,'sheet1','A1')
        xlswrite(xlsname,hrowy,'sheet2','A1')
    end
    mt = {fname};
    ntx2 = [Burstiness IntraBurstFrequency.mean IntraBurstFrequency.sd...
        IntraBurstSpikeNumber.mean IntraBurstSpikeNumber.sd...
        BurstLength.mean BurstLength.sd BurstFrequency];
    nty2 = [IntraBurstFrequency.all' IntraBurstSpikeNumber.all' BurstLength.all'];
    ntz2 = [IsiMatrix.mean IsiMatrix.sd IsiMatrix.num];
    hrowz = cell(1,1+2*size(IsiMatrix.mean,2)+size(IsiMatrix.num,2));
    hrowz{2} = 'mean';
    hrowz{1+size(IsiMatrix.mean,2)+1} = 'sd';
    hrowz{1+2*size(IsiMatrix.mean,2)+1} = 'num';
    str = ['A' num2str(size(ntx,1)+2)];
    xlswrite(xlsname,mt,'sheet1',str)
    str = ['B' num2str(size(ntx,1)+2)];
    xlswrite(xlsname,ntx2,'sheet1',str)
    str = ['A' num2str(size(nty,1)+2)];
    xlswrite(xlsname,mt,'sheet2',str)
    str = ['B' num2str(size(nty,1)+2)];
    xlswrite(xlsname,nty2,'sheet2',str)
    str = ['A' num2str(size(atz,1)+1)];
    xlswrite(xlsname,hrowz,'sheet3',str)
    str = ['A' num2str(size(atz,1)+2)];
    xlswrite(xlsname,mt,'sheet3',str)
    str = ['B' num2str(size(atz,1)+2)];
    xlswrite(xlsname,ntz2,'sheet3',str)
end
cd(mm)



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist2(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-11) '.mat'];
    end
end
files2 = files2(2:end);