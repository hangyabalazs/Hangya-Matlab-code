%% 

links2 = links;


%%

figure
[H,T]=dendrogram(links,0);
find(T==29)


%%

dlt = [109         776        1383        1474];
inx = [];
alldlt = [];
m = size(links,1);
while ~isempty(dlt)
    x = dlt(1);
    ins = find(links(:,1)==x);
    dlt = [dlt m+ins];
    inx = [inx ins];
    ins = find(links(:,2)==x);
    dlt = [dlt m+ins];
    inx = [inx ins];
    alldlt = [alldlt x];
    dlt = unique(setdiff(dlt,alldlt));
end
links(:,inx) = [];
links(inx,:) = [];

%%

K = flipud(str2num(get(gca,'YTickLabel')));
R = nan(30,2);
cntr = 0;
for k = K'
    cntr = cntr + 1;
    R(cntr,:) = [unique(c(find(T==k))) length(c(find(T==k)))];
end
R(6,:) = [];
