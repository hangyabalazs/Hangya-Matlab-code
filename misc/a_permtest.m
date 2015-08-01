function p = a_permtest(gene_no,Ano,Bno,olp)
%A_PERMTEST   Permutation test for overlap of genes.

% Permutation test, bootstrap sampling
bno = 1000;
genes = 1:gene_no;
olps = nan(1,bno);
for k = 1:bno
    rp = randperm(gene_no);
    Ac = genes(rp(1:Ano));
    rp = randperm(gene_no);
    Bc = genes(rp(1:Bno));
    olps(k) = length(intersect(Ac,Bc));
end

% Plot distribution
figure
hist(olps,50)
disp([mean(olps) median(olps)])

% p value
p = sum(olps>olp) / bno;