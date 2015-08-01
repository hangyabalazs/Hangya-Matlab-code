%% theta-delta

% size(delta_theta_alpha_single_hilb)
% ans =
%    13     3     4     3   100
% dim1=13 azaz 13 ksz
% dim2=3 azaz 3 frekvencia, 1-delta, 2- theta, 3- alpha
% dim3=4 azaz 4 koncicio, 1 -10%, 2 -37%, 3 -64% 4 -91%
% dim4=3 azaz 3 eeg csatorna, 1-Fz, 2-Cz, 3-Pz
% dim5=100, azaz 100 trial

data = delta_theta_alpha_single_hilb;
data(data<-pi) = data(data<-pi) + 2 * pi;
data(data>pi) = data(data>pi) - 2 * pi;
theta_fz_91 = cell(1,13);
edges = -pi:2*pi/9:pi;
cnts = (edges(1:end-1) + edges(2:end)) / 2;
figure
hold on
for k = 1:13
    theta_fz_91{k} = squeeze(data(k,2,4,1,:));
    [nm xout] = histc(theta_fz_91{k},edges);
    nm = nm(1:end-1);
    plot([cnts cnts+2*pi],[nm' nm'],'r')
end
delta_fz_91 = cell(1,13);
for k = 1:13
    delta_fz_91{k} = squeeze(data(k,1,4,1,:));
    [nm xout] = histc(delta_fz_91{k},edges);
    nm = nm(1:end-1);
    plot([cnts cnts+2*pi],[nm' nm'],'b')
end

%% theta-delta 2

edges2 = -pi:2*pi/18:pi;
cnts2 = (edges2(1:end-1) + edges2(2:end)) / 2;
figure
hold on
theta_fz_91_all = squeeze(data(:,2,4,1,:));
theta_fz_91_all = theta_fz_91_all(:);
[nm xout] = histc(theta_fz_91_all,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm' nm'],'r')
delta_fz_91_all = squeeze(data(:,1,4,1,:));
delta_fz_91_all = delta_fz_91_all(:);
[nm xout] = histc(delta_fz_91_all,edges2);
nm = nm(1:end-1);
plot([cnts2 cnts2+2*pi],[nm' nm'],'b')

[U2 p] = b_watsontwo(theta_fz_91_all,delta_fz_91_all)
b_circular_mean3(theta_fz_91_all) * 180 / pi
b_circular_mean3(delta_fz_91_all) * 180 / pi

%% theta-delta 3

p = [];
for k1 = 1:13
    for k2 = k1+1:13
        [U2 P] = b_watsontwo(theta_fz_91{k1},theta_fz_91{k2});
        p(end+1) = P(1);
    end
end
length(find(p<0.001))/length(p)