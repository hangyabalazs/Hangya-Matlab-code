%MEMORY_TEST    Test file for MATLAB memory usage.

next = 1;
while next < 39
    str = ['w' num2str(next) ' = [ones(10000,1000)];'];
    eval(str);
    next = next + 1
end
s = 0;