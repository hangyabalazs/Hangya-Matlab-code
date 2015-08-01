function [nb]= bayescount(totaltrials,binprob);
%compute the number of bins for bias calculation, using the
%Bayes counting procedure of Panzeri&Treves 96
%Arguments: 
%totaltrials: the number of trials in total
%binprob: the vector(1,totalbins) with the occupancy probabilities 
%totalbins: the total number of possible response bins is worked out from the size of the binprob argument

       
%gg is gamma
% nb is R_s
% xtr is iterated in order to compute the estimate of the expected
% no. of bins, and it becomes the number of "extra" relevant bins,
% Rs^tilde -Rs at the end. 
% qc_x is the Bayes estimate of the probabilities
% nbx is the expected no. of occupied bins

totalbins = length(binprob);   

nb = sum(binprob>eps);

if (nb < totalbins)
        nb_x = nb - sum((binprob>eps).*(binprob<1).*exp(log(1-binprob+eps)*totaltrials));
    delta_N_prev = totalbins;
    delta_N = abs(nb - nb_x);
    xtr = 0;
    while (delta_N<delta_N_prev & ((nb+xtr)<totalbins)) 
        xtr=xtr+1;
        nb_x = 0.0;
        gg = xtr*(1.-((totaltrials/(totaltrials+nb))^(1./totaltrials)));
        qc_x = (1-gg) * (binprob*totaltrials+1) / (totaltrials+nb);
        nb_x = sum((binprob>eps).*(1-exp(log(1-qc_x)*totaltrials)));
        qc_x = gg / xtr;
        nb_x = nb_x + xtr*(1. - exp(log(1. - qc_x) * totaltrials));
        delta_N_prev = delta_N;
        delta_N = abs(nb - nb_x);
    end; %while
    nb = nb + xtr - 1;
    if (delta_N < delta_N_prev)
        nb=nb+1;
    end   
end; %if (nb < totalbins)
