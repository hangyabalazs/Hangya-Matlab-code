function imicontrol1
%IMICONTROL1   Arteficial control for shifted mutual information calculation.
%   IMICONTROL1 generates arteficial control data by shifting a gaussian
%   window. See the comments in the code for the precise structure of the
%   data generated!
%
%   See also IMISHIFT.

lw = 800;
gs = gausswin(lw) * 10;           % 0.8 s gaussian 'wave'
data = zeros(20,10005);           % 10 second arteficial data, sampled on 1000 Hz
data(1,2001:2000+lw) = gs;        % wave at 2 sec. on channel1
data(2,2801:2800+lw) = gs;        % ch1 -> ch2 & ch6 (800 ms)
data(6,3101:3100+lw) = gs;
data(7,3701:3700+lw) = gs;        % ch2 & ch6 -> ch7 (900 ms)
data(8,4551:4550+lw) = gs;        % ch7 -> ch8 (850 ms)
data(12,4851:4850+lw) = gs;       % ch7 -> ch12 (850 ms)
data(20,5301:5300+lw) = gs;       % ch12 -> ch20 (750 ms)
data(19,6101:6100+lw) = gs;       % ch20 -> ch19 (800 ms)
data(15,6801:6800+lw) = gs;       % ch19 -> ch15 (700 ms)
data(20,7501:7500+lw) = gs;       % ch15 -> ch20 (700 ms)
data(10,8101:8100+lw) = gs;       % ch20 -> ch10 (600 ms)
data(20,8701:8700+lw) = gs;       % ch10 -> ch20 (600 ms)
data = data + rand(size(data))/5;   % 2% adding noise
data = data';
save('X:\In_Vivo\raw_data\human_SO\oiti37_lukacs\grid\mat_EEG_12\arteficial_control2.mat','data') % save