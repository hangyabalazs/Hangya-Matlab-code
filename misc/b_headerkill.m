function b_headerkill(filename)
%HEADERKILL Clears header from raw data file
%   HEADERKILL clears all nonnumeric lines from data file but do not separate the channels.
%
%   See also TEXTREAD.

% Measure runtime
tic

% Input argument check
error(nargchk(1,1,nargin));

% Change directory
cd d:\ms_neuron_rawdata\

% Clear headers
fid1 = fopen(filename,'r');
[pth fnm ext] = fileparts(filename);
new = [fnm '_V2.mat'];
str = ['fid2 = fopen(''' new ''',''w'');'];
eval(str)
next = 0;
c = 0;
op{1} = [];
while ~feof(fid1)
    next = next + 1;
    line = fgetl(fid1);
%     result1 = ~(isempty(findstr(line,'CHANNEL')));
%     result2 = ~(isempty(findstr(line,'volt')));
%     result3 = ~(isempty(findstr(line,'START')));
%     if result1 | result2 | result3
%         disp(next)
%     end
    if isempty(find(isletter(line)))
        if next ~= 1
            fprintf(fid2,'\n');
        end
        line = str2num(line);
        fprintf(fid2,'%f',line);
        line = str2num(line);
        op{c} = [op{c}; line];
    else
        c = c + 1;
        op{c} = [];
    end
end
fclose(fid1);
fclose(fid2);

% len = zeros(1,c);
% for i = 1:c
%     len(i) = length(op{i});
% end
% mlen = max(len);
% data = zeros(mlen,c);
% for i = 1:c
%     for j = 1:len(i)
%         data(j,i) = op{i}(j);
%     end
% end

% str = ['save ' new ' data'];
% eval(str);

% Measure runtime
toc