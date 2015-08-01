function intan_sync2(TTLs,saved_ts,x,y)

% load('C:\Balazs\courses\2013_TENSS\animal3\group_1_3_freely_moving1_processed2\position_timestamps.mat')

% intan_sync(TTLs,testxy_gr13(:,3),testxy_gr13(:,1),testxy_gr13(:,2))

son2 = TTLs;
ts = saved_ts / 1000;

% son2=son2(118:end);
% ts=ts(95:end);
% x=x(95:end);
% y=y(95:end);
% son2=son2(108:end);
% ts=ts(107:end);
% x=x(107:end);
% y=y(107:end);

% maxlen = min([length(son2) length(ts)]);
% son2 = son2(1:maxlen);
% ts = ts(1:maxlen);
% x = x(1:maxlen);
% y = y(1:maxlen);

t=ts-ts(1)+son2(1);


while any(abs(son2-(ts-ts(1)+son2(1)))>0.015)
    badinx = find(abs(son2-(ts-ts(1)+son2(1)))>0.015,1,'first');
    if son2(badinx)-(ts(badinx)-ts(1)+son2(1)) < 0.015
        son2(badinx) = [];
        ts(end) = [];
        x(end) = [];
        y(end) = [];
    else
        ts(badinx) = [];
        x(badinx) = [];
        y(badinx) = [];
        son2(end) = [];
    end
    
end

keyboard





% 
% % Match timestamps - in case of mismatch, try to fix
% 
%     % note: obsolete due the introduction of TTL parsing
% %     son2 = clearttls(son2); % eliminate recorded TTL's within 0.5s from each other - broken TTL pulse
%     if ~ismatch(ts,son2)
%         son2 = trytomatch(ts,son2);  % try to match time series by shifting
%         if ~ismatch(ts,son2)
%             son2 = tryinterp(ts,son2); % interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
%             if ~ismatch(ts,son2)  % TTL matching failure
%                 error('MakeTrialEvents:TTLmatch','Matching TTLs failed.')
%             else
%                 warning('MakeTrialEvents:TTLmatch','Missing TTL interpolated.')
%             end
%         else
%             warning('MakeTrialEvents:TTLmatch','Shifted TTL series.')
%         end
%     else
%         warning('MakeTrialEvents:TTLmatch','Broken TTLs cleared.')
%     end
% 
% % Eliminate last TTL's recorded in only one system
% sto = TE2.StimulusOn;
% if length(son2) > length(ts)   % time not saved in behavior file (likely reason: autosave was used)
%     son2 = son2(1:length(ts));
% elseif length(son2) < length(ts)  % time not recorded on Neuralynx (likely reason: recording stopped)
%     shinx = 1:length(son2);
%     ts = ts(shinx);
%     sto = sto(shinx);
%     TE2 = shortenTE(TE2,shinx);
% end
% TE2.TrialStart = son2 - sto;
% 
% % Save synchronized 'TrialEvents' file
% if ~isempty(TE2.TrialStart),
%     save([sessionpath filesep 'TrialEvents.mat'],'-struct','TE2')
% else
%     error('MakeTrialEvents:noOutput','Synchronization process failed.');
% end
% 
% % -------------------------------------------------------------------------
% function I = ismatch(ts,son2)
% 
% % Check if the two time series match notwithstanding a constant drift
% clen = min(length(ts),length(son2));
% I = abs(max(diff(ts(1:clen)-son2(1:clen)))) < 0.001;  % the difference between the timestamps on 2 systems may have a constant drift, but it's derivative should still be ~0
% 
% % note: abs o max is OK, the derivative is usually a small neg. number due
% % to drift of the timestamps; max o abs would require a higher tolerance
% % taking the drift into account (if 2 event time stamps are far, the drift
% % between them can be large)
% 
% % -------------------------------------------------------------------------
% function son2 = tryinterp(ts,son2)
% 
% % Interpolate missing TTL's or delete superfluous TTL's up to 10 erroneous TTl's
% for k = 1:10
%     if ~ismatch(ts,son2)
%         son3 = son2 - son2(1) + ts(1);
%         adt = diff(ts(1:min(length(ts),length(son2)))-son2(1:min(length(ts),length(son2))));
%         badinx = find(abs(adt)>0.1,1,'first') + 1;  % find problematic index
%         if adt(badinx-1) < 0    % interploate
%             ins = ts(badinx) - linterp([ts(badinx-1) ts(badinx+1)],[ts(badinx-1)-son3(badinx-1) ts(badinx+1)-son3(badinx)],ts(badinx));
%             son2 = [son2(1:badinx-1) ins+son2(1)-ts(1) son2(badinx:end)];
%         else
% %             ins = son3(badinx) - linterp([son3(badinx-1) son3(badinx+1)],[son3(badinx-1)-ts(badinx-1) son3(badinx+1)-ts(badinx)],son3(badinx));
% %             ts = [ts(1:badinx-1) ins ts(badinx:end)];
%             son2(badinx) = [];   % delete
%         end
%     end
% end
% 
% % -------------------------------------------------------------------------
% function son2 = trytomatch(ts,son2)
% 
% % Try to match time series by shifting
% len = length(son2) - 15;
% minx = nan(1,len);
% for k = 1:len
%     minx(k) = max(diff(ts(1:15)-son2(k:k+14)));  % calculate difference in the function of shift
% end
% mn = min(abs(minx));
% minx2 = find(abs(minx)==mn);
% minx2 = minx2(1);   % find minimal difference = optimal shift
% son2 = son2(minx2:min(minx2+length(ts)-1,length(son2)));
% 
% % -------------------------------------------------------------------------
% function son2 = clearttls(son2)
% 
% % Eliminate recorded TTL's within 0.5s from each other
% inx = [];
% for k = 1:length(son2)-1
%     s1 = son2(k);
%     s2 = son2(k+1);
%     if s2 - s1 < 0.5
%         inx = [inx k+1]; %#ok<AGROW>
%     end
% end
% son2(inx) = [];
% 
% % -------------------------------------------------------------------------
% function TE2 = shortenTE(TE2,shinx)
% 
% % Eliminate behavioral trials
% fnm = fieldnames(TE2);
% for k = 1:length(fieldnames(TE2))
%     TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
% end