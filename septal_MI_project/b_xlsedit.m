function b_xlsedit
%XLSEDIT   Reads xls files and unites content in proper order.
%   XLSEDIT copies data from an xls file to another using first column to
%   specify order. Edit code to fit XLSEDIT to specific problem!
%
%   See also XLSREAD.

% File names
global DATAPATH
fn1 = [DATAPATH 'Entry3b\Stat\power\noth_rorc.xls'];
fn2 = [DATAPATH '2Uentropy\3summary16b.xls'];

% Read
sh1 = 'NothNoburst';
sh2 = 'power_ol1';
[N1, T1, rawdata1] = xlsread(fn1,sh1,'A2:M60');
[N2, T2, rawdata2] = xlsread(fn2,sh2,'A4:AA86');

% Write
T22 = T2(:,1);
for i = 1:length(T1)
    name = T1{i};
    inx = find(strcmp(name(1:6),T22));
    for j = 1:12
        rawdata2{inx,15+j} = rawdata1{i,1+j};
    end
end

xlswrite(fn2,rawdata2,sh2,'A4');