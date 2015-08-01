%SGEPOCHPHASE2   Creates delta phase data.
%   SGEPOCHPHASE2 calculates delta phase values for filtered and epoched
%   data from Experiment2 by the use of the Hilbert-transform.
%
%   See also SGFILT5CHANS_FILTART2, SGEPOCH2 and SGPCONTROLPHASE2.

for k1 = 1:11
    for k2 = 1:size(erpdata21{k1},3)
        ed = erpdata21{k1}(1,:,k2);
        eds = angle(hilbert(ed));
        sh21FZ(k1).d(k2) = eds(925);
        ed = erpdata21{k1}(2,:,k2);
        eds = angle(hilbert(ed));
        sh21CZ(k1).d(k2) = eds(925);
        ed = erpdata21{k1}(3,:,k2);
        eds = angle(hilbert(ed));
        sh21PZ(k1).d(k2) = eds(925);
    end
    for k2 = 1:size(erpdata41{k1},3)
        ed = erpdata41{k1}(1,:,k2);
        eds = angle(hilbert(ed));
        sh41FZ(k1).d(k2) = eds(925);
        ed = erpdata41{k1}(2,:,k2);
        eds = angle(hilbert(ed));
        sh41CZ(k1).d(k2) = eds(925);
        ed = erpdata41{k1}(3,:,k2);
        eds = angle(hilbert(ed));
        sh41PZ(k1).d(k2) = eds(925);
    end
end

save X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2filtart/filtartdelta.mat ...
    sh21FZ sh21CZ sh21PZ sh41FZ sh41CZ sh41PZ