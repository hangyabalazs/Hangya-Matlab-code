%AGANNOT   Annotation class for ANALYSE_SNIFFWHISK GUI.
%   AGANNOT defines classe for user annotations. 
%
%   See also ANALYSE_SNIFFWHISK.

classdef agannot < handle
    properties (SetAccess = private)
        StartPoint = NaN;
        EndPoint = NaN;
        AnnotString = '';
        LogDirectory
        IsSaved = 0;
    end
    properties (SetAccess = private, Transient)
        StartLineHandle = struct;
        EndLineHandle = struct;
        PatchHandle = struct;
        SerialNumber
        ParrentGui
        ParrentAxes = struct;
        IsDeleted = 0;
    end
    methods
        
        % Constructor
        function obj = agannot(sp,ep,str,handles)
            obj.StartPoint = sp;
            obj.EndPoint = ep;
            obj.AnnotString = str;
            obj.SerialNumber = length(handles.annot) + 1;
            gui_pos = get(handles.figure1,'Position');
            G = agannot_gui(gui_pos,handles.figure1);   % evoke GUI for interactive setting
            uiwait(G)
            global ANNOTOUT
            annotout = ANNOTOUT;
            obj.StartPoint = annotout.sp;
            obj.EndPoint = annotout.ep;
            obj.AnnotString = annotout.str;
            obj.ParrentGui = handles.figure1;
            global DATAPATH
            obj.LogDirectory = [DATAPATH 'analyse_gui\annotations\' handles.filename];
            if ~isdir(obj.LogDirectory)
                mkdir(obj.LogDirectory)
            end
            
            fid = fopen([obj.LogDirectory '\log.txt'],'a');   % write log file
            fprintf(fid,'%f %f %f %f %f %f \n',clock);
            fprintf(fid,'%f %f %s %f %f \n',obj.StartPoint,obj.EndPoint,...
                obj.AnnotString,obj.SerialNumber,obj.IsDeleted);
            fclose(fid);
        end
        
        % Overload isempty for annot objects
        function I = isempty(obj)
            if (isempty(obj.StartPoint) || isnan(obj.StartPoint)) && ...
                    (isempty(obj.EndPoint) || isnan(obj.EndPoint)) && ...
                    isempty(obj.AnnotString)
                I = 1;
            else
                I = 0;
            end
        end
        
        % Overload plot for annot objects
        function plot(obj,y_lim,str)
            obj.ParrentAxes.(str) = get(obj.ParrentGui,'CurrentAxes');
            hcmenu = uicontextmenu;
            hcedit = 'editannot(getappdata(gco,''annotobj''))';
            hcdelete = 'delline(getappdata(gco,''annotobj''))';
            uimenu(hcmenu,'Label','edit','Callback',hcedit);
            uimenu(hcmenu,'Label','delete','Callback',hcdelete);
            if ~isempty(obj.StartPoint)     % plot a line for start point
                L = line([obj.StartPoint obj.StartPoint],y_lim,'Color',...
                    [0.1 0.9 0.1],'LineWidth',2);
                setappdata(L,'annotobj',obj)
                obj.StartLineHandle.(str) = L;
                set(L,'uicontextmenu',hcmenu)
            end
            if ~isempty(obj.EndPoint)     % plot a line for end point
                L = line([obj.EndPoint obj.EndPoint],y_lim,'Color',...
                    [0.25 0.75 0.25],'LineWidth',2);
                setappdata(L,'annotobj',obj)
                obj.EndLineHandle.(str) = L;
                set(L,'uicontextmenu',hcmenu)
            end
            if ~isempty(obj.StartPoint) && ~isempty(obj.EndPoint)     % plot a patch for interval
                P = patch([obj.EndPoint obj.EndPoint obj.StartPoint obj.StartPoint obj.EndPoint],...
                    [y_lim y_lim(2) y_lim(1) y_lim(1)],[0 1 0],'FaceColor',[0.25 0.75 0.25],'FaceAlpha',0.5);
                setappdata(P,'annotobj',obj)
                obj.PatchHandle.(str) = P;
                set(P,'uicontextmenu',hcmenu)
            end
        end
        
        % Set parrent GUI for loaded objects
        function setparrentgui(obj,parrentgui)
            obj.ParrentGui = parrentgui;
            obj.StartLineHandle = struct;
            obj.EndLineHandle = struct;
            obj.PatchHandle = struct;
            obj.ParrentAxes = struct;
        end
        
        % Edit an existing annotation
        function editannot(obj)
            gui_pos = get(obj.ParrentGui,'Position');
            G = agannot_gui(gui_pos,obj.ParrentGui,obj.StartPoint,...
                obj.EndPoint,obj.AnnotString);   % evoke GUI for interactive setting
            uiwait(G)
            global ANNOTOUT
            annotout = ANNOTOUT;
            if (isempty(annotout.sp) || isnan(annotout.sp)) && ...
                    (isempty(annotout.ep) || isnan(annotout.ep)) && ...
                    isempty(annotout.str)     % return if editing is cancelled
                return
            end
            obj.StartPoint = annotout.sp;
            obj.EndPoint = annotout.ep;
            obj.AnnotString = annotout.str;
            obj.IsSaved = 0;
            str = fieldnames(obj.StartLineHandle);
            y_lim = {};
            bdf = {};
            for k = 1:length(str)
                y_lim{k} = get(obj.StartLineHandle.(str{k}),'YData');
                bdf{k} = get(obj.StartLineHandle.(str{k}),'ButtonDownFcn');   % get buttondown function (zoom)
            end
            
            delline(obj)   % delete from GUI
            for k = 1:length(str)
                set(obj.ParrentGui,'CurrentAxes',obj.ParrentAxes.(str{k}));
                plot(obj,y_lim{k},str{k})   % replot changed annotation
                set(obj.StartLineHandle.(str{k}),'ButtonDownFcn',bdf{k});   % assign buttondown function (zoom)
                set(obj.EndLineHandle.(str{k}),'ButtonDownFcn',bdf{k});
                set(obj.PatchHandle.(str{k}),'ButtonDownFcn',bdf{k});
            end
            obj.IsDeleted = 0;
            
            fid = fopen([obj.LogDirectory '\log.txt'],'a');   % write log file
            fprintf(fid,'%f %f %f %f %f %f \n',clock);
            fprintf(fid,'%f %f %s %f %f \n',obj.StartPoint,obj.EndPoint,...
                obj.AnnotString,obj.SerialNumber,obj.IsDeleted);
            fclose(fid);
        end
        
        % Remove from plots and mark as deleted
        function delline(obj)
            hans = cell2mat(struct2cell(obj.StartLineHandle));
            delete(hans(ishandle(hans)))
            hans = cell2mat(struct2cell(obj.EndLineHandle));
            delete(hans(ishandle(hans)))
            hans = cell2mat(struct2cell(obj.PatchHandle));
            delete(hans(ishandle(hans)))
            obj.StartLineHandle = struct;
            obj.EndLineHandle = struct;
            obj.PatchHandle = struct;
            obj.IsDeleted = 1;
            
            fid = fopen([obj.LogDirectory '\log.txt'],'a');   % write log file
            fprintf(fid,'%f %f %f %f %f %f \n',clock);
            fprintf(fid,'%f %f %s %f %f \n',obj.StartPoint,obj.EndPoint,...
                obj.AnnotString,obj.SerialNumber,obj.IsDeleted);
            fclose(fid);
        end
        
        % Delete
        function delete(obj)
            disp('deleting')
            hans = cell2mat(struct2cell(obj.StartLineHandle));
            delete(hans(ishandle(hans)))
            hans = cell2mat(struct2cell(obj.EndLineHandle));
            delete(hans(ishandle(hans)))
            hans = cell2mat(struct2cell(obj.PatchHandle));
            delete(hans(ishandle(hans)))
        end
        
        % Save
        function obj = saveobj(obj)
            obj.IsSaved = 1;    % set 'IsSaved' property
        end
    end
end