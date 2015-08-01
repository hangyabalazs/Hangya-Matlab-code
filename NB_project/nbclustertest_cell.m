%%

L = [];
CN = [];
C = {};
for k = 2:8
    cmbs = nchoosek(1:8,k);
    
    for t = 1:size(cmbs,1);
        dmtxstr = [];
        for st = 1:k
            dmtxstr = [dmtxstr 'v' num2str(cmbs(t,st)) ' '];
        end
        
        eval(['dmtx = [' dmtxstr '];']);
        dist = pdist(dmtx);
        links = linkage(dist,'Ward');
        c = cluster(links,4);
        pChATvinx = setdiff(find(c==1),ChATvinx);   % putative cholinergic neurons, new indexing
        inxset = [ChATactinx; NTactinx];   % indices for index conversion
        pChATactinx = inxset(pChATvinx);   % putative cholinergic neurons, old indexing
        pChAT = tags(pChATactinx);   % cell IDs of pChAT cells
        NTvinx2 = setdiff(NTvinx,pChATvinx);   % remove pChAT from NT
        C1inx = sort(inxset(c==1));
        C2inx = sort(inxset(c==2));
        C3inx = sort(inxset(c==3));
        C4inx = sort(inxset(c==4));
        chatno1 = length(intersect(ChATactinx,C1inx));
        chatno2 = length(intersect(ChATactinx,C2inx));
        chatno3 = length(intersect(ChATactinx,C3inx));
        chatno4 = length(intersect(ChATactinx,C4inx));
        if chatno1 + chatno2 == 0 |...
                chatno1 + chatno3 == 0 |...
                chatno1 + chatno4 == 0 |...
                chatno2 + chatno3 == 0 |...
                chatno2 + chatno4 == 0 |...
                chatno3 + chatno4 == 0
            [jnk chatcluster] = max([chatno1 chatno2 chatno3 chatno4]);
            alll = [length(C1inx) length(C2inx) length(C3inx) length(C4inx)];
            L(end+1) = alll(chatcluster);
            C{end+1} = cmbs(t,:);
            CN(end+1) = jnk;
        end
    end
end