function b = adjust_track(a,end_position) 

b = (a - a(1)) / (a(end) - a(1)) * (end_position - a(1)) + a(1);

% c = b;
% c(2:end) = sqrt(c(2:end).^2-(1580-1460).^2);
% 
% d = (c - c(1)) / (c(adjusted_index) - c(1)) * (end_position - c(1)) + c(1);
% 
% e = (c - c(1)) / (c(adjusted_index) - c(1)) * ((end_position - offset) - c(1)) + c(1) + offset;