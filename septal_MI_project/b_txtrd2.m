function [fname,f1,f2] = b_txtrd(t);
%TXTRD   Gets the data out of an "_ivph_result" type text file.  
%   The input argument is the "serial number" of the file you want to examine.
%
%   See also TEXTREAD.

%Input arguments check
error(nargchk(1,1,nargin));

%Textread
t = ['f:\_analysis\analyse_result\ivph_texts\_ivph_result' num2str(t) '.txt'];
fname = textread(t,'%s');
fname = fname(1);
[f1, f2] = textread(t,'%f%f','headerlines',1);