function [tI cI] = b_wread2
%WREAD2   Import Tisza data.
%
%   See also WDISC and IN4.

% Read
fname = 'e:\Tisza\Mellekfolyok\tbecsq.txt';
tbc = wrd(fname);
tbc = dfill(tbc)

% -------------------------------------------------------------------------
function data = wrd(fname)

[data1 data2 data3 data4 data5 data6] = textread(fname,'%f %f %f %f %f %f');
data = [data1 data2 data3 data4 data5 data6];

fnd1 = find(data1==2000,1,'first');
fnd2 = find(data2(fnd1:end)==7,1,'first');
i_first = fnd1 + fnd2 - 1;
fnd1 = find(data1==2002,1,'first');
fnd2 = find(data2(fnd1:end)==12,1,'first');
fnd3 = find(data3(fnd1+fnd2-1:end)==30,1,'first');
i_second = fnd1 + fnd2 + fnd3 - 2;
data(i_first,:)
data(i_second,:)

data;

% -------------------------------------------------------------------------
function data = dfill(data)

data1 = data(:,1);
data2 = data(:,2);
data3 = data(:,3);
data4 = data(:,4);
data5 = data(:,5);
data6 = data(:,6);
data1_new = [];
data2_new = [];
data3_new = [];
data4_new = [];
data5_new = [];
data6_new = [];

df = diff(data5);
fn = find(df~=15);
for k = length(data1)
    d =  df(k);
    if isequal(mod(d,60),15)
        data1_new(end+1) = data1(k);
        
        if d < 0
            h = 1;
        else
            h = 0;
        end
        if ~isequal(data4(k)-data4(k-1),h)
            error(['More than an hour missing at index' num2str(k) '.'])
        end
        
    elseif isequal(mod(d,60),30)
        data1_new(end+1) = (data1(k-1) + data1(k)) / 2;
        data1_new(end+1) = data1(k);
        
        if d < 0
            h = 1;
        else
            h = 0;
        end
        if ~isequal(data4(k)-data4(k-1),h)
            error(['More than an hour missing at index' num2str(k) '.'])
        end
        
    elseif isequal(mod(d,60),45)
        data1_new(end+1) = (2 * data1(k-1) + data1(k)) / 3;
        data1_new(end+1) = (data1(k-1) + 2 * data1(k)) / 3;
        data1_new(end+1) = data1(k);
        
        if d < 0
            h = 1;
        else
            h = 0;
        end
        if ~isequal(data4(k)-data4(k-1),h)
            error(['More than an hour missing at index' num2str(k) '.'])
        end
                
    elseif isequal(mod(d,60),0)
        data1_new(end+1) = (3 * data1(k-1) + data1(k)) / 4;
        data1_new(end+1) = (data1(k-1) + data1(k)) / 2;
        data1_new(end+1) = (data1(k-1) + 3 * data1(k)) / 4;
        data1_new(end+1) = data1(k);
        
        if ~isequal(data4(k)-data4(k-1),1)
            error(['More than an hour missing at index' num2str(k) '.'])
        end
        
    else
        error(['Time difference between measure poits:' num2str(d) ' at index' num2str(k) '.'])
    end
end
    
