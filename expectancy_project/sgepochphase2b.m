%SGEPOCHPHASE2B   Creates delta phase data.
%   SGEPOCHPHASE2B calculates delta phase values for unfiltered and epoched
%   data from Experiment2 by the use of the Hilbert-transform (after
%   filtering).
%
%   See also SGFILT5CHANS_FILTART2, SGEPOCH2B and SGPCONTROLPHASE2.

flt = fir1(256,[30 70]/250); %4096
for k1 = 1:11
    for k2 = 1:size(erpdata21{k1},3)
        ed = erpdata21{k1}(1,:,k2);
        ed = filtfilt(flt,1,ed);
        eds = angle(hilbert(ed));
        sh21FZ(k1).d(k2) = eds(925);
        ed = erpdata21{k1}(2,:,k2);
        ed = filtfilt(flt,1,ed);
        eds = angle(hilbert(ed));
        sh21CZ(k1).d(k2) = eds(925);
        ed = erpdata21{k1}(3,:,k2);
        ed = filtfilt(flt,1,ed);
        eds = angle(hilbert(ed));
        sh21PZ(k1).d(k2) = eds(925);
    end
    for k2 = 1:size(erpdata41{k1},3)
        ed = erpdata41{k1}(1,:,k2);
        ed = filtfilt(flt,1,ed);
        eds = angle(hilbert(ed));
        sh41FZ(k1).d(k2) = eds(925);
        ed = erpdata41{k1}(2,:,k2);
        ed = filtfilt(flt,1,ed);
        eds = angle(hilbert(ed));
        sh41CZ(k1).d(k2) = eds(925);
        ed = erpdata41{k1}(3,:,k2);
        ed = filtfilt(flt,1,ed);
        eds = angle(hilbert(ed));
        sh41PZ(k1).d(k2) = eds(925);
    end
end

save X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2nofilt/filtartgamma.mat ...
    sh21FZ sh21CZ sh21PZ sh41FZ sh41CZ sh41PZ