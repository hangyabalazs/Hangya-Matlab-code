%%

ChAT = {'n071_141129a_2.3' 'n071_141201a_2.3' 'n071_141201a_7.2' 'n071_141218a_6.1'...
'n071_141229a_3.2' 'n072_141222a_4.3' 'n072_141223a_4.2'}
cellid = ChAT{1};

%%

cellids{1}='n067_141009b_2.2'
cellids{2}='n067_141011x_1.1'
cellids{3}='n067_141012a_3.2'
cellids{4}='n067_141017a_1.3'
cellids{5}='n067_141017a_2.1'
cellids{6}='n067_141017a_4.1'
cellids{7}='n067_141018a_5.1'
cellids{8}='n067_141019a_5.2'
cellids{9}='n070_141112a_4.1'
cellids{10}='n070_150104a_5.1'
cellids{11}='n070_150106a_5.2'
cellids{12}='n078_150104a_1.1'
cellids{13}='n078_150110a_3.1'

%%

cellids=ChAT;
for k=1:length(cellids)
    cellid=cellids{k}
    [psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
        'DeliverFeedback',[-0.6 0.6],...
        'dt',0.001,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput','FalseAlarm==1','maxtrialno',Inf,...
        'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
    log10(stats.maxvalue./stats.baseline)
end

%%

[psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
    'DeliverFeedback',[-1 1],...
    'dt',0.001,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
    'event_filter','custom','filterinput','FalseAlarm==1','maxtrialno',Inf,...
    'baselinewin',[-0.5 0],'testwin',[0 0.1],'relative_threshold',0.1);