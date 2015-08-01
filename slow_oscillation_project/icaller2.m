function icaller2
%ICALLER2   Caller function for SO project functions.
%
%   See also IMISIGNIO.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients_all.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call
for k = 11:26
    disp(k)
    pat = atx{k,1};
    patno = num2str(atx{k,2});
    eg = num2str(atx{k,3});
    nm_rows = atx{k,4};
    nm_cols = atx{k,5};
    esg = [atx{k,6} atx{k,7}];
    
%     ifftsnr2(pat,patno,eg,nm_rows,nm_cols,esg)
    imisignio_lin2_stand(pat,patno,eg,nm_rows,nm_cols,esg)
%     imisignio_stand(pat,patno,eg,nm_rows,nm_cols,esg)
    imisignio_stand(pat,patno,eg,nm_rows,nm_cols,esg)
    
    close all
end