%%

sp=speed;
al=arrowlength_orig;
as=arrowstrength;

%%

sp(2,:)=speed;
al(2,:)=arrowlength_orig;
as(2,:)=arrowstrength;

%%

sp(3,:)=speed;
al(3,:)=arrowlength_orig;
as(3,:)=arrowstrength;

%%

figure
bar(sp,'stacked')
axis off

%%

figure
hold on
plot(asx,as(1,:)/sum(as(1,:)),'b')
plot(asx,as(2,:)/sum(as(2,:)),'r')
plot(asx,as(3,:)/sum(as(3,:)),'m')

%%

figure
hold on
stairs(alx,al(1,:),'b')
stairs(alx,al(2,:),'r')
stairs(alx,al(3,:),'m')