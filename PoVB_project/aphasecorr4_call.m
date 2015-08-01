function aphasecorr4_call
%APHASECORR4_CALL    Calls APHASE_CORR4B for a sequence of directories.
%
%   See also APHASECORR_CALL and APHASE_CORR4B.

% Directories
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call APHASE_CORR4B
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    aphase_corr4b(inpdir)
end