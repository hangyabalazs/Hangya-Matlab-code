%NLXEVENTS   Event class for ANALYSE_SNIFFWHISK.
%   NLXEVENTS defines classes for Neuralynx events. It plots lines for
%   events.
%
%   See also ANALYSE_SNIFFWHISK.

classdef NlxEvents < handle
    properties (SetAccess = private)
        EventID
        TimeStamp
        Nttl
        EventString
    end
    properties (SetAccess = private, Transient)
        LineHandle = struct;
        TextHandle1 = [];
        TextHandle2 = [];
    end
    methods
        function obj = NlxEvents(eid,ts,nttl,estr)
            if nargin > 0
                obj.EventID = eid;
                obj.TimeStamp = ts;
                obj.Nttl = nttl;
                obj.EventString = estr;
            end
        end
        function plot(obj,y_lim,str)
            L = line([obj.TimeStamp obj.TimeStamp],y_lim,'Color','r');
            obj.LineHandle.(str) = L;
        end
        function delline(obj)
            hans = cell2mat(struct2cell(obj.LineHandle));
            delete(hans(ishandle(hans)))
            hans = cell2mat(struct2cell(obj.TextHandle1));
            delete(hans(ishandle(hans)))
            hans = cell2mat(struct2cell(obj.TextHandle2));
            delete(hans(ishandle(hans)))
            obj.LineHandle = struct;
            obj.TextHandle1 = struct;
            obj.TextHandle2 = struct;
        end
        function delete(obj)
            hans = cell2mat(struct2cell(obj.LineHandle));
            delete(hans(ishandle(hans)))
        end
    end
end