%SGEPOCH2   Epoching (from S.G.).
%   SGEPOCH2 creates epoched data from filtered raw data from Experiment2.
%
%   See also SGFILT5CHANS_FILTART2, SGEPOCHPHASE2 and SGCONTROLPHASE2.

filename=dir('X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2filtart/*.set');
filename=struct2cell(filename); filename=filename(1,:);
triggers = { '21', '41'};
for f=1:size(filename,2)

EEG=pop_loadset(strcat('X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2filtart/',char(filename(f)))); % cd('D:\ATTENTION_HIGHGAMMA\setdataconcat');D:\ATTENTION_HIGHGAMMA\setdataconcatFZCZPZ
%EEG = pop_select( EEG, 'channel',[2 4 5 6 8] ); %EEG = pop_select( EEG,'channel',[2:3:8] );
EEG = eeg_checkset( EEG ); % 2:3:


EEGep = pop_epoch( EEG, triggers(1), [-0.5  4.05],   'epochinfo','yes');
erpdata21{f} = EEGep.data; clear EEGep;
EEGep = pop_epoch( EEG, triggers(2), [-0.5  4.05],   'epochinfo','yes');
erpdata41{f} = EEGep.data; clear EEGep;

end
save X:\In_Vivo\raw_data\human_SG\expectancy_rawDATA_EXP2filtart/filtartdata.mat erpdata21 erpdata41