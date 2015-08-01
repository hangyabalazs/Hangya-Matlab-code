%SGEPOCH   Epoching (from S.G.).
%   SGEPOCH creates epoched data from filtered or unfiltered raw data.
%
%   See also SGFILT5CHANS_FILTART, SGEPOCHPHASE and SGPHASE_FILTART3.


filename=dir('X:\In_Vivo\raw_data\human_SG\setdataconcatFZCZPZC3C4filtart_ms1500/*.set');
filename=struct2cell(filename); filename=filename(1,:);
triggers = { '102', '372', '642', '912' };
for f=1:size(filename,2)

EEG=pop_loadset(strcat('X:\In_Vivo\raw_data\human_SG\setdataconcatFZCZPZC3C4filtart_ms1500/',char(filename(f)))); % cd('D:\ATTENTION_HIGHGAMMA\setdataconcat');D:\ATTENTION_HIGHGAMMA\setdataconcatFZCZPZ
%EEG = pop_select( EEG, 'channel',[2 4 5 6 8] ); %EEG = pop_select( EEG,'channel',[2:3:8] );
EEG = eeg_checkset( EEG ); % 2:3:


for tr = 1:4
     EEGep = pop_epoch( EEG, triggers(tr), [-1.5  1.5],   'epochinfo','yes');
     if size(EEGep.data,3) > 100
         EEGep.data = EEGep.data(:,:,1:100);
     end
     erpdata(f,tr,:,:,:) = EEGep.data; clear EEGep;
end


end
save X:\In_Vivo\raw_data\human_SG\setdataconcatFZCZPZC3C4filtart_ms1500/erpdata.mat erpdata