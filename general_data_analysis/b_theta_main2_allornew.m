function b_theta_main2
%THETA_MAIN2    Ultimate theta selector.
%   THETA_MAIN2 calls the thetaselector functions. It performs the entire theta
%   segment selection.
%
%   THETA_MAIN2 works on two directories: the raw data directory and the results'
%   directory. You are able to modify them through editing the program code.
%
%   See also THETASELECTORRUN2, THETA and NONTHETA.

% Input argument check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH     % b_startup is required!
global DATADIR
inpdir = [DATADIR 'raphe_matfiles\temp\discriminated\'];
resdir = [DATAPATH 'Wavelet_raphe\'];
cd(resdir)

% All or new
aon = input('All files in directory or new ones only? [a/n] ','s');
if isequal(aon,'n')
    k = dir(inpdir);
    lenlist = length(k)-2;
    liststring = {};
    for t = 1:lenlist
        if ~k(t+2).isdir
            liststring{end+1} = k(t+2).name(1:6);
        end
    end
    if ~b_isdir2('theta_segments')
        readystring = {};
    else
        k2 = dir([resdir 'theta_segments\']);
        lenlist2 = length(k2)-2;
        readystring = {};
        for t = 1:lenlist2
            if ~k2(t+2).isdir
                readystring{end+1} = k2(t+2).name(16:21);
            end
        end
    end
    list = setdiff(liststring,readystring);
    cd(inpdir)
    mkdir('temp')       % create temporary folder
    for z = 1:length(list)
        str=['copyfile(''' list{z} '*'',''temp'',''f'');'];
        eval(str);
    end
    inpdir = fullfile(inpdir,'temp\');
    cd(resdir)
    mkdir('temp')       % create temporary folder
    resdir = fullfile(resdir,'temp\');
elseif isequal(aon,'a')
else
    error('Unexpected input.')
end

% Thetaselectorrun
disp('Running THETASELECTORRUN2...');
b_thetaselectorrun2(inpdir,resdir);

% Load scalevector and newstep
load scalevector.mat
load newstep.mat

% Theta
disp('Running THETA...');
b_theta(f,newstep,resdir)

% Nontheta
disp('Running NONTHETA...')
b_nontheta(f,newstep,resdir,inpdir)

% Remove temporary folder
if isequal(aon,'n')
    rmdir(inpdir,'s')
end