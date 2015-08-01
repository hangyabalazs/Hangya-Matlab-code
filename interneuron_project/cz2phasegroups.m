%%

[ntz mtz atz] = xlsread('X:\In_Vivo\balazs\_analysis\Czurko2\VMfit\EM_grouping.xls','sheet2');


%%

inx1 = find(ntz(:,1)==1);
phase1 = ntz(inx1,2);
phase1_rad = phase1 / 180 * pi;
min(phase1)
max(phase1)
min(phase1_rad)
max(phase1_rad)
circular_mean(phase1,'deg')

%%

inx2 = find(ntz(:,1)==2);
phase2 = ntz(inx2,2);
phase2_rad = phase2 / 180 * pi;
min(phase2)
max(phase2)
min(phase2_rad)
max(phase2_rad)
circular_mean(phase2,'deg')

%%

inx3 = find(ntz(:,1)==3);
phase3 = ntz(inx3,2);
phase3_rad = phase3 / 180 * pi;
min(phase3)
max(phase3)
min(phase3_rad)
max(phase3_rad)
circular_mean(phase3,'deg')

%%

inx4 = find(ntz(:,1)==4);
phase4 = ntz(inx4,2);
phase4_rad = phase4 / 180 * pi;
min(phase4)
max(phase4)
min(phase4_rad)
max(phase4_rad)
circular_mean(phase4,'deg')

%%

phase_rad = ntz(:,2) / 180 * pi;
X1 = vmpdf(phase_rad,mu(1),kappa(1));
X2 = vmpdf(phase_rad,mu(2),kappa(2));
X3 = vmpdf(phase_rad,mu(3),kappa(3));
X4 = vmpdf(phase_rad,mu(4),kappa(4));

X = [X1 X2 X3 X4];
for k = 1:74
    mG(k)=find(X(k,:)==max(X(k,:)));
end
mG2 = mG';

%%

[p,table,stats,terms] = anovan(ntz(:,5),ntz(:,1))

[p,table,stats,terms] = anovan(ntz(:,8),ntz(:,1))
[p,table,stats,terms] = anovan(ntz(:,9),ntz(:,1))
[p,table,stats,terms] = anovan(ntz(:,10),ntz(:,1))

[p,table,stats,terms] = anovan(ntz(:,11),ntz(:,1))

%%

length(find(ntz(:,1)==1&ntz(:,7)==1))
length(find(ntz(:,1)==1&ntz(:,7)==0))

length(find(ntz(:,1)==2&ntz(:,7)==1))
length(find(ntz(:,1)==2&ntz(:,7)==0))

length(find(ntz(:,1)==3&ntz(:,7)==1))
length(find(ntz(:,1)==3&ntz(:,7)==0))

length(find(ntz(:,1)==4&ntz(:,7)==1))
length(find(ntz(:,1)==4&ntz(:,7)==0))

