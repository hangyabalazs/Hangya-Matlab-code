%%

bst = length(fun_nostim);
for k = 1:bst
    pk_nostim(k) = fun_nostim{k}.mu;
    pk_stim(k) = fun_stim{k}.mu;
    pk_diff(k) = fun_nostim{k}.mu - fun_stim{k}.mu;
    
    wdt_nostim(k) = fun_nostim{k}.sigma;
    wdt_stim(k) = fun_stim{k}.sigma;
    
    sh_nostim(k) = fun_nostim{k}.b;
    sh_stim(k) = fun_stim{k}.b;
    
    sc_nostim(k) = fun_nostim{k}.a;
    sc_stim(k) = fun_stim{k}.a;
end

%%

sum(pk_stim>fun_nostim_orig.mu) / bst
sum(wdt_stim>fun_nostim_orig.sigma) / bst
sum(sh_stim>fun_nostim_orig.b) / bst
sum(sc_stim>fun_nostim_orig.a) / bst