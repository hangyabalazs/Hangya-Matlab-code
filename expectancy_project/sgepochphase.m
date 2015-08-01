%SGEPOCHPHASE   Creates delta phase data.
%   SGEPOCHPHASE calculates delta phase values for filtered and epoched
%   data by the use of the Hilbert-transform.
%
%   See also SGFILT5CHANS_FILTART, SGEPOCH and SGPHASE_FILTART3.

flt = fir1(256,[0.5 3]/250); %4096
erpphase = zeros(13,4,5,1,100);
for k1 = 1:13
    for k2 = 1:4
        for k3 = 1:5
            for k5 = 1:100
                ed = erpdata(k1,k2,k3,:,k5);
%                 ed = filtfilt(flt,1,ed);
%                 ed(1:1499) = ed(1500);
                ah = angle(hilbert(ed));
                erpphase(k1,k2,k3,1,k5) = ah(1500);
            end
        end
    end
end

save X:\In_Vivo\raw_data\human_SG\setdataconcatFZCZPZC3C4filtart_ms1500/filtartdelta.mat erpphase