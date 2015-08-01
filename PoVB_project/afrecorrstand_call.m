function afrecorrstand_call
%AFRECORRSTAND_CALL    Calls AFRE_CORR_STAND for a sequence of directories.
%
%   See also APHASECORR_CALL and AFRE_CORR_STAND.

% Directories
inproot = 'Y:\extra32\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call AFRE_CORR_STAND
for k = 1:length(inpadd)
    cn = str2num(inpadd{k}(4:5));
    if cn > 30
        inpdir = [inproot inpadd{k} '\bas\']
        afre_corr_stand(inpdir)
    end
end