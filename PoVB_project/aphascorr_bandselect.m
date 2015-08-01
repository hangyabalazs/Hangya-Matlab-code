function aphasecorr_bandselect
%APHASECORR_BANDSELECT   EEG frequency - phase angle correlation on pooled data.
%   APHASECORR_BANDSELECT loads EEG frequency - mean_angle data and creates
%   a cumulative plot for frequency band selection.
%
%   See also APHASE_CORR4B.

% Directories
global DATAPATH
inpdir = [DATAPATH 'Andi\Ketxyl\PhaseCorr4b_3\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\FreBandSelect\'];
mm = pwd;

% Filelist
[files files_short] = filelist2(inpdir);
sf = length(files_short);

% Main
Ang_bas = [];
Ang_bic = [];
Fre_bas = [];
Fre_bic = [];
for o = 1:sf
    fname = files(o).name;     % load
    if ~isempty(findstr(fname,'ANGFRECORR.mat'))
        load([inpdir fname])
        Ang_bas = [Ang_bas meanang_bas];
        Ang_bic = [Ang_bic meanang_bic];
        Fre_bas = [Fre_bas fre_bas];
        Fre_bic = [Fre_bic fre_bic];
    end
end

% Plot and save
cd(resdir)
H = figure;
plot(Ang_bas,Fre_bas,'.')   % mean angle vs. EEG frequency
hold on
plot(Ang_bic,Fre_bic,'r.')
xlabel('mean angle')
ylabel('EEG frequency')
saveas(H,'AngFreCorr_new.fig')
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