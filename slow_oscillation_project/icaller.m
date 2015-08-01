function icaller
%ICALLER   Caller function for SO project functions.
%
%   See also ICONTRIBMAP and IMISORUN_RANDCONVDIV.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call
for k = 1:15
    k
    pat = atx{k,1};
    patno = num2str(atx{k,2});
    eg = num2str(atx{k,3});
    nm_rows = atx{k,4};
    nm_cols = atx{k,5};
    
    imisorun_randconvdiv(pat,patno,eg,nm_rows,nm_cols)
    
    close all
end