%% F-test to compare SD values (expectancy project)

for k=1:13
    s10{k}=squeeze(FiveChannel3_phase_rt(k,1,2,3,:));
    s37{k}=squeeze(FiveChannel3_phase_rt(k,2,2,3,:));
    s64{k}=squeeze(FiveChannel3_phase_rt(k,3,2,3,:));
    s91{k}=squeeze(FiveChannel3_phase_rt(k,4,2,3,:));
    p10(k,1:100)=squeeze(FiveChannel3_phase_rt(k,1,2,3,:));
    p37(k,1:100)=squeeze(FiveChannel3_phase_rt(k,2,2,3,:));
    p64(k,1:100)=squeeze(FiveChannel3_phase_rt(k,3,2,3,:));
    p91(k,1:100)=squeeze(FiveChannel3_phase_rt(k,4,2,3,:));
end
figure;[nm xout]=hist(p10(:),15);stairs(xout,nm)
hold on;[nm xout]=hist(p37(:),15);stairs(xout,nm,'Color','g')
hold on;[nm xout]=hist(p64(:),15);stairs(xout,nm,'Color','c')
hold on;[nm xout]=hist(p91(:),15);stairs(xout,nm,'Color','r')

%%

for k=1:13
    sd10(k) = std(p10(k,:));
    sd37(k) = std(p37(k,:));
    sd64(k) = std(p64(k,:));
    sd91(k) = std(p91(k,:));
end

figure
plot(1,sd10,'b.')
hold on
plot(2,sd37,'b.')
plot(3,sd64,'b.')
plot(4,sd91,'b.')
xlim([0 5])

%%

for k=1:13
    se10(k) = std(p10(k,:)) / 10;
    se37(k) = std(p37(k,:)) / 10;
    se64(k) = std(p64(k,:)) / 10;
    se91(k) = std(p91(k,:)) / 10;
end

figure
plot(1,se10,'b.')
hold on
plot(2,se37,'b.')
plot(3,se64,'b.')
plot(4,se91,'b.')
xlim([0 5])

%%

[stat,h,p] = b_Ftest(sd10,sd37,0.05)
[stat,h,p] = b_Ftest(sd37,sd64,0.05)
[stat,h,p] = b_Ftest(sd64,sd91,0.05)

%%

for k = 1:13
    [stat,h,p] = b_Ftest(p10(k,:),p91(k,:),0.05)
end
