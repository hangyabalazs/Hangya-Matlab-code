function alignevent = findAlignEvent_negfeedback_gonogo(cellid)
%FINDALIGNEVENT_NEGFEEDBACK_GONOGO   Dynamic event defintion.
%   ALIGNEVENT = FINDALIGNEVENT_NEGFEEDBACK_GONOGO(CELLID) allows dynamic
%   defintion of trial events. The name of the negactive feedback event
%   (ALIGNEVENT) is returned for CELLID ('DeliverFeedback' or 'LeftPortIn'
%   for false alarms, depending on the session type).
%
%   See also FINDALIGNEVENT_POSFEEDBACK_GONOGO and ULTIMATE_PSTH.

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
    alignevent = 'LeftPortIn';
end