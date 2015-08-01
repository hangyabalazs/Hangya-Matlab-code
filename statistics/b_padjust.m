function p = b_padjust(p,method,n)

switch method
    case 'holm'
        i = (1:n);
        p = min(1,max((n - i + 1) * p));
    case 'hochberg'
        i = (n:-1:1);
        p = min(1,min((n - i + 1) * p));
    case 'fdr'
        i = (n:-1:1);
        p = min(1,min(n./i * p));
    case 'bonferroni'
        p = min(n*p,1);
end