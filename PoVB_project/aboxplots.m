%% boxplot

H = figure;
boxplot([b1' b2' b3' b4' b5'],[zeros(size(b1')) ones(size(b2')) 2*ones(size(b3')) ...
    3*ones(size(b4')) 4*ones(size(b5'))],'labels',[{'so_sw vs so_fw'} ...
    {'so_sw vs fo_sw'} {'so_sw vs fo_fw'} {'so_fw vs fo_fw'} {'fo_fw vs fo_sw'}]);
% [Wp_eu,Wh_eu] = b_ranksum2(m1,m2,'alpha',0.05);
% if Wh_eu
%     clr = 'red';
% else
%     clr = 'black';
% end
% y_lim = ylim;
% x_lim = xlim;
% tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
% tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 3 / 5;
% text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
% tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
% tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
% text(tpos1,tpos2,loc)

%% dots

hold on
plot(ones(size(b1)),b1,'.')
plot(2*ones(size(b2)),b2,'.')
plot(3*ones(size(b3)),b3,'.')
plot(4*ones(size(b4)),b4,'.')
plot(5*ones(size(b5)),b5,'.')

%% circ mean + circ SE

ldf = 2;
mndf = zeros(1,ldf);
SE = zeros(1,ldf);
for k = 1:ldf
    [ftm, ang, mvl] = mvlmn(B{k},'deg');
    mndf(k) = ang;
    SE(k) = circular_SE(B{k},'deg');
end
H = figure;
bar((1:ldf),mndf)
hold on
errorbar((1:ldf),mndf,SE,'k+')

%% circ mean

H = figure;
bar((1:ldf),mndf,'FaceColor','none')
hold on
for k = 1:ldf
    plot(k*ones(size(B{k})),B{k},'.')
end

%% boxplot

H = figure;
boxplot([b1' b2'],[zeros(size(b1')) ones(size(b2'))],'labels',[{'so_sw vs so_fw'} ...
    {'fo_fw vs fo_sw'}]);
% [Wp_eu,Wh_eu] = b_ranksum2(m1,m2,'alpha',0.05);
% if Wh_eu
%     clr = 'red';
% else
%     clr = 'black';
% end
% y_lim = ylim;
% x_lim = xlim;
% tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
% tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 3 / 5;
% text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
% tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
% tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
% text(tpos1,tpos2,loc)

%% dots

hold on
plot(ones(size(b1)),b1,'.')
plot(2*ones(size(b2)),b2,'.')

%% dots

hold on
plot(ones(size(b1)),a1(:,1),'.')
plot(2*ones(size(b2)),a2(:,2),'.')
