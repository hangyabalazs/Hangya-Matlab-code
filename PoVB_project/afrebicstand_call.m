function afrebicstand_call
%AFREBICSTAND_CALL    Calls AFE_BIC_STAND for a sequence of directories.
%
%   See also AFRE_BIC_STAND.

% Directories
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Read Excel file
fn = [tabledir 'tablazat_balazsnak.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call AFRE_BIC_STAND
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    afre_bic_stand2(inpdir,atx)
end