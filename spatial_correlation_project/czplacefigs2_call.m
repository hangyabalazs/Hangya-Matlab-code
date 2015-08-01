function czplacefigs2_call
%CZPLACEFIGS2_CALL   Caller function for CZPLACEFIGS2.
%   CZPLACEFIGS2_CALL calls CZPLACEFIGS2 for negative or positive pairs.
%
%   See also CZPLACEFIGS2.

% Directories
global DATAPATH
inpdir_xls = [DATAPATH 'Czurko\discriminated2\'];

% Load
type = 'neg';
xlsname = [inpdir_xls 'placefigs.xls'];
headerrows = 1;
[ntz mtz atz] = xlsread(xlsname,type);
mtz(1:headerrows,:) = [];
atz(1:headerrows,:) = [];
sf = size(mtz,1);   % number of pairs
for k = 1:sf
    [C Cb Cc trsr] = czplacefigs3(atz{k,2},atz{k,3},atz{k,1},type);
    op(k,1:3) = atz(k,1:3);
    op{k,4} = C;
    op{k,5} = Cb;
    op{k,6} = Cc;
    op{k,7} = trsr;
end

% Write excel file output
xlsname = [inpdir_xls 'placedata.xls'];
xlswrite(xlsname,op,type,'A1')