%CONSERVE   Contains useful parts for memory conserving programming technique.

% ------------------------------------------------------------------------
% Saving to temporary directory
% ------------------------------------------------------------------------

% Make sure tempdir is on the path
added_temp = checktempdir;

% Save variables in a temporary directory
tmpname = [tempname '.mat'];
%save 

% Be sure to clean up temp file
if exist(tmpname)
    delete(tmpname); 
end
if added_temp
    rmpath(tempdir)
end

% ------------------------------------------------------------------------
function added_temp = checktempdir

if isempty(findstr(tempdir,matlabpath))
  addpath(tempdir)
  added_temp = 1;
else
  added_temp = 0;
end
% Generate a temp m-file name that is used to hold evaluated strings
% for compilation.  Test temp file to verify it can be opened in write
% mode and deleted.  Unfortunately the delete command does not provide
% an error return value so I'm catching the warning and then erroring out.
lastwarn('');
tmpname = [tempname '.m'];
fid = fopen(tmpname, 'w');
fclose(fid);
if strcmp(lastwarn,'File not found or permission denied')
    error(['Unable to run WAVELET due to inability to write to and remove ' ...
            'temporary file (' tmpname ').']);
end
if exist(tmpname), 
    delete(tmpname);
end

% ------------------------------------------------------------------------
% Calculating power and phase
% ------------------------------------------------------------------------
sw1 = size(wave,1);
sw2 = size(wave,2);
pieceno = 10;
segm = fix(sw2/pieceno);
inx2 = 0;
next = 1;
while inx2 < sw2
    inx1 = max(1,inx2+1);
    inx2 = min(sw2,inx1+segm);
    wavefrag = wave(:,inx1:inx2);
    str = ['save temp' num2str(next) ' wavefrag'];
    eval(str)
    clear wavefrag
    next = next + 1;
end
clear wave
power = [];
phase = [];
for wsn = 1:next-1
    str = ['load temp' num2str(wsn)];
    eval(str)
    powerfrag = (abs(wavefrag)) .^2;
    power = [power powerfrag];
    clear powerfrag
    phasefrag = angle(wavefrag);
    phase = [phase phasefrag];
    clear phasefrag
    clear wavefrag
end
% ------------------------------------------------------------------------

sw2 = size(wave,2);
pieceno = 20;
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