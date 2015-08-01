function alignevent = findAlignEvent_posfeedback_gonogo(cellid)
%FINDALIGNEVENT_POSFEEDBACK_GONOGO   Dynamic event defintion.
%   ALIGNEVENT = FINDALIGNEVENT_POSFEEDBACK_GONOGO(CELLID) allows dynamic
%   defintion of trial events. The name of the positive feedback event
%   (ALIGNEVENT) is returned for CELLID ('DeliverFeedback' or
%   'LeftWaterValveOn' for hits, depending on the session type).
%
%   See also FINDALIGNEVENT_NEGFEEDBACK_GONOGO and ULTIMATE_PSTH.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   07-May-2012

%   Edit log: BH 5/20/13

% Checking whether 'DeliverFeedback' event is available
sesstype = getvalue('session_type',cellid);
if isequal(sesstype,{'feedbackdelay'})
    alignevent = 'DeliverFeedback';
else
    alignevent = 'LeftWaterValveOn';
end