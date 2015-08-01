function [tI cI] = b_wread
%WREAD   Import Tisza data.
%
%   See also WDISC and IN4.

% Read 'tivInnov'
tI1 = textread('e:\Tisza\tivInnov.txt','%s');
lentI = length(tI1);
tI = zeros(1,lentI);
for i = 1:lentI
    if isequal(tI1{i},'NA')
        if isequal(i,1)
            tI(i) = str2num(tI1{2});
        else
            tI(i) = (str2num(tI1{i-1}) + str2num(tI1{i+1})) / 2;
        end
    else
        tI(i) = str2num(tI1{i});
    end
end
tI = tI(2:end);     % tI is shifted!

% Read 'csenInnov'
cI1 = textread('e:\Tisza\csenInnov.txt','%s');
lencI = length(cI1);
cI = zeros(1,lencI);
for i = 1:lencI
    if isequal(cI1{i},'NA')
        if isequal(i,2)
            cI(i) = str2num(cI1{2});
        else
            cI(i) = (str2num(cI1{i-1}) + str2num(cI1{i+1})) / 2;
        end
    else
        cI(i) = str2num(cI1{i});
    end
end
cI = cI(1:end-1);