% get the list of sessions to convert
% convert them
% save the TE file
% tic;
basedir='/Users/ranades/SoloData/Data/hanna/';
animal_name={'hr020'};
% animal_name={'hr020','hr021','hr022'};
% animal_name={'hr013'};
animal_name={'a2'};
% % animal_name='a1';
% animal_name='a3';
% % 
basedir='/Users/ranades/SoloData/Data/jordan/';
% % % animal_name='jw02';
% animal_name={'jw01','jw02','jw03'};
% animal_name={'jw04','jw05','jw06','jw07','jw08'};
animal_name={'jw09','jw010','jw011','jw012','jw013','jw014','jw015','jw016','jw017','jw018','jw019'};
% animal_name={'jw05'};
% % % animal_name='a2';
% basedir='/Users/ranades/SoloData/Data/matt/';
% % % % animal_name='jw02';
% % animal_name={'jw01','jw02','jw03'};
% animal_name={'mk01','mk02','mk03','mk04','mk05'};
% animal_name='a1';
% animal_name='a3';
tic;
for iA=1:length(animal_name),
    allfiles=listfiles([basedir animal_name{iA}]);
    sessions=setdiff(allfiles,allfiles(strmatch('TE',allfiles)));
    for iS=1:length(sessions),
        filepath=[basedir animal_name{iA} filesep sessions{iS}];
        try
            TE=solo2trialevents2_auditory_gonogo(filepath);
            savename = ['TE_' animal_name{iA} '_' TE.sessionID{end}];
            save([basedir animal_name{iA} filesep savename],'-struct','TE')
        catch
        end
    end
end
toc;