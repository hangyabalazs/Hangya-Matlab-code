function [f,newstep] = b_thetaselectorrun_long(inpdir,resdir)
%THETASELECTORRUN_LONG   Version of THETASELECTORRUN for long registrations.
%   You have to use THETASELECTORRUN_LONG for registrations longer than 200 seconds.
%   THETASELECTORRUN saves a binary file called 'long.mat', which contains the identifiers
%   of the long registrations. You have to put the "long files" in the input directory
%   of THETASELECTORRUN_LONG (either manually or using LONGCOPY). For details, see 
%   THETASELECTORRUN.
%
%   THETASELECTORRUN_LONG uses two directories: one for the data files and one for the 
%   results. You are able to modify these directories through editing the program code,
%   or you can add them as input parameters: THETASELECTORRUN_LONG(INPUT_DIRECTORY,
%   RESULT_DIRECTORY). Input directory can be specified using THETASELECTORRUN_LONG
%   (INPUT_DIRECTORY) syntax.
%
%   [F,NEWSTEP] = THETASELECTORRUN_LONG(INPDIR,RESDIR) returns F scalevector and NEWSTEP
%   downsampling constant (step for downsamlinp 10 kHz sampled EEG).
%
%   See also THETASELECTOR3, THETASELECTOR_BETA3, LONGCOPY and THETASELECTORRUN.

% Input arguments check
error(nargchk(0,2,nargin));

% Define directories
global DATAPATH
if isempty(DATAPATH)
    clear glogal DATAPATH;
    b_startup
    global DATAPATH
end
global DATADIR
if isempty(DATADIR)
    clear glogal DATADIR;
    b_startup
    global DATADIR
end
if nargin < 2
    resdir = [DATAPATH, 'Wavelet2\'];  %Here are the results
end
if nargin < 1
    inpdir = [DATAPATH,'DATA\analysenow4\'];    %Here are the data files
end
cd(resdir)
create_subdir       %create subdirectories

% Import
files = dir(inpdir);
files = files(3:end);
sf = length(files);
mm = pwd;

wb = waitbar(0,'Running THETASELECTORRUN_LONG...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

for o = 1:sf
    for dtx = 1:2
        fname = files(o).name;
        ffnm = [inpdir fname];
        data = b_load_data(ffnm);
        meret = size(data,1);
        
% Computing the input variables
        if dtx == 1
            datinx1 = 1; %first point of the interval
            datinx2 = 2000000;  %last point of the interval
        else
            datinx1 = 2000001;  %first point of the interval
            datinx2 = size(data,1); %last point of the interval
        end
        if datinx2 - datinx1 > 2000000
            global LONG
            LONG{end+1} = fname(1:6);
            continue
            waitbar(o/sf)   %Progress indicator
        end
        b_imp(fname,inpdir,data,datinx1,datinx2)
        
% Free memory
        save temp o sf wb files inpdir dtx
        clear
        clear function
        
% Theta selection
        [wave,f,newstep] = b_waveletcall_for_applications;
        
        sw2 = size(wave,2);
        pieceno = 5;
        segm = fix(sw2/pieceno);
        power = [];
        while ~isempty(wave)
            index1 = 1;
            index2 = min(segm,size(wave,2));
            wavefrag = wave(:,index1:index2);
            powerfrag = (abs(wavefrag)) .^ 2;
            clear wavefrag
            wave(:,index1:index2) = [];
            power = [power powerfrag];
        end
        clear wave
        
        [H,OM] = b_thetaselector3(power,f,newstep);
        ratio = b_thetaselector_beta3(power,f,newstep);
        
% Saving
        Out = zeros(3,sw2);
        Out(1:2,:) = OM;
        Out(3,:) = ratio;
        
        global IN
        fname = IN{3};
        pont = findstr(fname,'.');
        filenam = fname(1:pont(1)-1);
        fig = gcf;
        ax = findobj(fig,'Type','axes');
        axes(ax(2))
        title(['EEG ',filenam(1:3),' ',filenam(5:6)]);
        load temp
        if dtx == 1
            cd jpg
            eval(['saveas(H,''THETA_SELECT_LONG1_',filenam(1:6),'.jpg'')']);
            cd ..
            cd matrix
            eval(['save(''THETA_SELECT_LONG1_',filenam(1:6),'.mat'',''Out'')']);
            cd ..
        else
            cd jpg
            eval(['saveas(H,''THETA_SELECT_LONG2_',filenam(1:6),'.jpg'')']);
            cd ..
            cd matrix
            eval(['save(''THETA_SELECT_LONG2_',filenam(1:6),'.mat'',''Out'')']);
            cd ..
        end
    end
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator
delete temp.mat
close all

% Save 'long'
global LONG
long = LONG;
save long long
clear global LONG

% Save scalevector and 'newstep'
save scalevector f
save newstep newstep

% ----------------------------------------------------------------------------------
function create_subdir

% Create subdirectories
if ~b_isdir2(thetaselection_long)
    mkdir(thetaselection_long)
end
cd thetaselection_long
if ~b_isdir2('matrix')
    mkdir matrix
end
if ~b_isdir2('jpg')
    mkdir jpg
end