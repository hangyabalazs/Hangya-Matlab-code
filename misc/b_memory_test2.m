% function b_memory_test
%MEMORY_TEST    Test file for MATLAB memory usage.

% z = [];
% next = 0;
% while next < 1000
%     z = [z rand(1,1000)];
%     next = next + 1;
% end

% z = zeros(100,100,100,30);
% assignin('base','MemTestOutput',z);

next = 1;
while 1
    str = ['w' num2str(next) ' = [ones(10000,1000)];'];
    eval(str);
    next = next + 1
end
s = 0;