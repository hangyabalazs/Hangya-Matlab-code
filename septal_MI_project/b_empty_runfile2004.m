function OutputArg = b_empty_runfile2004(InputArg)
% EMPTY_RUNFILE2004   Example RunFile for EMPTY_GUI2004.
%   Input argument is a cell array of the input arguments specified in EMPTY_GUI2004.
%   All output should be stored in one output argument.
%
%   See also EMPTY2004 and EMPTY_GUI2004.

x = InputArg{1};
y = sin(x);
H = figure;
plot(x,y);
OutputArg = struct('y',y,'H',H);

global FID
fid = FID;
fprintf(fid,'%f',H);