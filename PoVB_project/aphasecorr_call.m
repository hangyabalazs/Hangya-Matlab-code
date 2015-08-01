function aphasecorr_call
%APHASECORR_CALL    Calls APHASE_CORR2B for a sequence of directories.
%
%   See also APHASE_CORR2, APHASE_CORR3 and APHASE_CORR2B.

% Directories
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call APHASE_CORR2B
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\bic\']
    aphase_corr2b(inpdir)
end    