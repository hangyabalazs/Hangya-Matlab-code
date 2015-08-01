%%

tt = 8;
% for tt=69:99
fld = fieldnames(segments3);
tmpp = segments3(tt);
for k = 8:11
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
tmpp.end = tmpp.exhalation_end(end);
tmp = tmpp.(fld{12});
tmp(tmp<tmpp.start|tmp>tmpp.end) = [];
tmpp.(fld{12}) = tmp;
segments3(tt) = tmpp;
% end



%% sort

lense = length(segments3);
fr = nan(1,lense);
for k = 1:lense
   fr(k) = segments3(k).inhalation_start(2) - segments3(k).start;
end
[sr ins] = sort(fr);
segments3 = segments3(ins);
segments3 = fliplr(segments3);

%%

tt = 10;
% for tt=1:58
fld = fieldnames(segments3);
tmpp = segments3(tt);
tmpp.start = tmpp.inhalation_start(1);
tmpp.end = tmpp.exhalation_end(end);
tmp = tmpp.(fld{12});
tmp(tmp<tmpp.start|tmp>tmpp.end) = [];
tmpp.(fld{12}) = tmp;
segments3(tt) = tmpp;
% end