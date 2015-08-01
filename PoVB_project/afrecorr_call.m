function afrecorr_call
%AFRECORR_CALL    Calls AFRE_CORR for a sequence of directories.
%
%   See also APHASECORR_CALL and AFRE_CORR.

% Directories
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call AFRE_CORR
for k = 1:length(inpadd)
    cn = str2num(inpadd{k}(4:5));
    if cn > 30
        inpdir = [inproot inpadd{k} '\bas\']
        afre_corr(inpdir)
    end
end