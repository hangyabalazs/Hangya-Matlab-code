%% sgrtstat

[s1 s2] = size(allrt);
mn = mean(allrt);
SE = std(allrt) ./ sqrt(s1);
figure
hold on
bar(1:s2,mn)
errorbar(1:s2,mn,SE,'k+')
ylim([220 320])
set(gca,'YTick',[220:20:320])
set(gca,'XTick',[1:4])
set(gca,'XTickLabel',{'10%' '37%' '64%' '91%'})
set(gca,'FontSize',16)
xlim([0.5 4.5])

%% boxplot

figure
boxplot(allrt)

%% Wilcowon signrank-test

[Wp_eu,Wh_eu] = b_signrank2(allrt(:,1),allrt(:,2),'alpha',0.05)
[Wp_eu,Wh_eu] = b_signrank2(allrt(:,2),allrt(:,3),'alpha',0.05)
[Wp_eu,Wh_eu] = b_signrank2(allrt(:,3),allrt(:,4),'alpha',0.05)