function [H,slope] = somslope2(DataPoints,lim1,lim2)

% Slice
n = length(DataPoints);



wn = DataPoints;
linlen = lim2 - lim1 + 1;
lsrate = 1000;
time = linspace(0,linlen/lsrate,linlen);
[gr,icp,err] = ...
        linefit(time*1000,wn(lim1:lim2));    % time is rescaled because large scaling difference screws up linefit
    gr = gr * 1000;    % rescale to agree with original time vector
% r(gr<0) = -1000;   % looking for rising edge, therefore set error of descending segments to -inf
slope = gr;
intercept = icp;

H = [];
H = figure;
plot(wn)
hold on
t = lim1:lim2;
y = time .* slope + intercept;
plot(t,y,'r')
err