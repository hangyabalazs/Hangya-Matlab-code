function saveas2(varargin)
%SAVEAS2 Saves figure with actual size.
%
%   SYNTAX:
%     saveas2
%     saveas2(NAME)                           % NAME may be [] 
%     saveas2(NAME,RESOLUTION)
%     saveas2(NAME,RESOLUTION,FORMAT)         % FORMAT may be: 'unknown'
%     saveas2(NAME,RESOLUTION,FORMAT,PRINTOPT)
%     saveas2(HANDLE,...)
%
%   INPUT:
%     HANDLE     - Figure handle. See NOTE below.
%                  DEFAULT: gcf
%     NAME       - String specifying the output figure name file. May
%                  include full/relative path and with/without extension.
%                  DEFAULT: Promts gui for user input if empty.
%     RESOLUTION - Figure resolution.
%                  DEFAULT: 150 (150 dpi).
%     FORMAT     - Output figure format. Use 'unknown' to display a list.
%                  See NOTE below. 
%                  DEFAULT '.png'
%     PRINTOPT   - Extra options for the PRINT function.
%                  DEFAULT: none
%
%    DESCRIPTION:
%      Depending on the FORMAT, this functions uses PRINT or SAVEAS
%      functions to save the specified (HANDLE) figure(s) matching the size
%      it has on screen.
%
%      This is a mayor modification of Gabe Hoffmann function SAVE2PDF that
%      allows multiple formats (including '.fig'). Besides, it may display
%      a gui for output filename and image extension.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * When specifing several handles, the same quantity of filenames
%       shold be given as a cell of strings.
%     * When no FORMAT is given, it is extracted from the NAME extension.
%     * If FORMAT was invalid, the user is asked to select a valid one from
%       a list. See PRINTTABLE for details.
%     * If the filename already exist it will be OVERWRITTEN without any
%       warning!!!
%
%   EXAMPLE:
%     figure, imagesc(peaks(40))
%     % Saves current figure asking for the name
%       saveas2
%     % Saves with specific name, format and resolution
%       saveas2('peaks.jpg',300)
%     % Saves the '.fig' file
%       saveas2('peaks.fig')
%     % Saves displaying a list of extensions (or devices)
%       saveas2('peaks',[],'unknown')
%    
%   SEE ALSO:
%     SAVEAS, PRINT
%     and
%     SAVE2PDF by Gabe Hoffmann
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   saveas2.m
%   VERSION: 1.1 (Sep 30, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Aug 20, 2009)
%   1.1      Fixed important bug with 'fig' format and multiple figures.
%            (Sep 30, 2009) 

%   DISCLAIMER:
%   saveas2.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Defaults:
NAME       = '';
HANDLE     = gcf;
FORMAT     = 'png';
RESOLUTION = '150';
PRINTOPT   = {};

% Parameters:
givenFORMAT = false;
dformat     = FORMAT;

% First checks HANDLE.
if ~isempty(varargin) && any(ishandle(varargin{1}))
 HANDLE = varargin{1};
 varargin(1) = [];
end

% Now other inputs.
if ~isempty(varargin)
 % Checks NAME.
 if ~isempty(varargin{1}) && ischar(varargin{1})
  NAME = varargin{1}; 
 end
 varargin(1) = [];
 if ~isempty(varargin)
  % Checks RESOLUTION.
  if ~isempty(varargin{1}) && isnumeric(varargin{1})
   RESOLUTION = int2str(varargin{1});
  end
  varargin(1) = [];
  if ~isempty(varargin)
   % Checks FORMAT.
   if ~isempty(varargin{1}) && ischar(varargin{1})
    FORMAT = varargin{1}(1,:);
    givenFORMAT = true;
    if strcmp(FORMAT(1),'.')
     FORMAT(1) = '';
     if isempty(FORMAT)
      givenFORMAT = false;
     end
    end    
   end
   varargin(1) = [];
   if ~isempty(varargin)
    % Checks PRINTOPT.
    PRINTOPT = varargin;
    clear varargin
   end
  end
 end
end

% Gets size.
N = length(HANDLE);

% Checks NAME:
if isempty(NAME)
 if N~=1
  error('CVARGAS:saveas2:tooManyHandles',...
   'NAMEs should be given when more than one figure.')
 end
 [thefile,thedir] = uiputfile('*','Save to file:');
 if thefile==0; return; end
 NAME = fullfile(thedir,thefile);
end

% Changes char to cell.
if ischar(NAME)
 NAME = {NAME};
end

% Checks sizes.
Nn = length(NAME);
if  Nn~=N
 error('CVARGAS:saveas2:incorrectHandlesAndNamesLengths',...
   'Length NAME and HANDLE must match.')
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------
 
% Main loop
checkDev = true;
theext = '';
for cont = 1:N
 
 name   = NAME{cont};
 handle = HANDLE(cont);
 if ~ishandle(handle)
  continue
 end
 
 % Looks for extension for FORMAT.
 if ~givenFORMAT
  [thedir,thename,theext] = fileparts(name);
  if length(theext)>2
   FORMAT = theext(2:end);
  else
   FORMAT = dformat;
  end
  name = fullfile(thedir,thename);
 end
 
 % Uses SAVEAS.
 if strcmpi(FORMAT,'m')    || strcmpi(FORMAT,'fig') || ...
    strcmpi(FORMAT,'mfig') || strcmpi(FORMAT,'mmat')
   saveas(handle,name,FORMAT)
   continue % Fixed BUG Sep 2009
 end
 
 % Backup original figure settings.
 prePaperType     = get(handle,'PaperType');
 prePaperUnits    = get(handle,'PaperUnits');
 preUnits         = get(handle,'Units');
 prePaperPosition = get(handle,'PaperPosition');
 prePaperSize     = get(handle,'PaperSize');
 
 % Make changing paper type possible
 set(handle,'PaperType','<custom>'); 
 % Set units to all be the same
 set(handle,'PaperUnits','inches');
 set(handle,'Units'     ,'inches');
 % Set the page size and position to match the figure's dimensions
 position = get(handle,'Position');
 set(handle,'PaperPosition',[0 0 position(3:4)]);
 set(handle,'PaperSize'    ,position(3:4));
 
 % Checks print device.
 if checkDev
  % Gets valid devices.
  [ops,dev,ext] = printtables;
  [ext,ind] = sort(ext);
  dev = dev(ind);
  % Looks extensions.
  i = strmatch(FORMAT,ext);
  if length(i) >= 1
   i = i(1);
  elseif isempty(i)
   % Looks devices.
   i = strmatch(FORMAT,dev,'exact');
   if isempty(i)
    i = strmatch(FORMAT,dev);
    if ~isempty(i)
     i = i(1);
    end
   end
  end
  if ~isempty(i)
   FORMAT = dev{i};
  else
   % Promts for print device.
   ind = menu('SELECT IMAGE FILE EXTENSION (or DEVICE)',...
      cellfun(@(a,b)[upper(a) '  (' b ')'],ext,dev,'UniformOutput',false));
   if ind==0
    return
   end
   FORMAT = dev{ind};
   name = [name theext];
   if ~isempty(ext{ind}) && ~strcmp(theext,ext{ind})
    name = [name '.' ext{ind}];
   end  
  end
  if ~givenFORMAT
   checkDev = false;
  end
 end
 
 % Uses PRINT.
 print(handle,name,['-d' FORMAT],['-r' RESOLUTION],PRINTOPT{:})
 
 % Restores figure settings.
 set(handle,'PaperType'    ,prePaperType);
 set(handle,'PaperUnits'   ,prePaperUnits);
 set(handle,'Units'        ,preUnits);
 set(handle,'PaperPosition',prePaperPosition);
 set(handle,'PaperSize'    ,prePaperSize);
  
end


% [EOF]   saveas2.m