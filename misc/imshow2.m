function h=imshow2(varargin)
%IMSHOW Display image.  
%   IMSHOW(I) displays the intensity image I.
%
%   IMSHOW(I,[LOW HIGH]) displays I as a grayscale intensity image,
%   specifying the display range for I. The value LOW (and any value less
%   than LOW) displays as black, the value HIGH (and any value greater than
%   HIGH) displays as white, and values in between display as intermediate
%   shades of gray. IMSHOW uses the default number of gray levels. If you
%   use an empty matrix ([]) for [LOW HIGH], IMSHOW uses [min(I(:))
%   max(I(:))]; the minimum value in I displays as black, and the maximum
%   value displays as white.
%
%   IMSHOW(RGB) displays the truecolor image RGB.
%
%   IMSHOW(BW) displays the binary image BW. Values of 0 display as black, and
%   values of 1 display as white.
%
%   IMSHOW(X,MAP) displays the indexed image X with the colormap MAP.
%
%   IMSHOW(FILENAME) displays the image contained in the graphics file
%   FILENAME.  The file must contain an image that can be read by IMREAD or
%   DICOMREAD. IMSHOW calls IMREAD or DICOMREAD to read the image from the
%   file, but the image data is not stored in the MATLAB workspace.  If the
%   file contains multiple images, the first one will be displayed. The file
%   must be in the current directory or on the MATLAB path.
%
%   HIMAGE = IMSHOW(...) returns the handle to the image object created by
%   IMSHOW.
%
%   IMSHOW(...,PARAM1,VAL1,PARAM2,VAL2,...)  specifies parameters and
%   corresponding values that control various aspects of the image display.
%   Parameter names can be abbreviated, and case does not matter.
%
%   Parameters include:
%
%   'DisplayRange'           Two-element vector [LOW HIGH] that controls the
%                            display range of an intensity image. See above
%                            for more details about how to set [LOW HIGH].
%
%                            Including the parameter name is optional, except
%                            when the image is specified by a filename. 
%                            The syntax IMSHOW(I,[LOW HIGH]) is equivalent to
%                            IMSHOW(I,'DisplayRange',[LOW HIGH]).
%                            The parameter name must be specified when 
%                            using IMSHOW with a filename, as in the syntax
%                            IMSHOW(FILENAME,'DisplayRange'[LOW HIGH]).
%
%   'InitialMagnification'   A numeric scalar value, or the text string 'fit',
%                            that specifies the initial magnification used to 
%                            display the image. When set to 100, the image is 
%                            displayed at 100% magnification. When set to 
%                            'fit' IMSHOW scales the entire image to fit in 
%                            the window.
%
%                            On initial display, the entire image is visible.
%                            If the magnification value would create an image 
%                            that is too large to display on the screen,  
%                            IMSHOW warns and displays the image at the 
%                            largest magnification that fits on the screen.
%
%                            By default, the initial magnification is set to
%                            the value returned by
%                            IPTGETPREF('ImshowInitialMagnification').
%
%                            If the image is displayed in a figure with its
%                             'WindowStyle' property set to 'docked', then
%                            IMSHOW warns and displays the image at the
%                            largest magnification that fits in the figure.
%
%   'XData'                  Two-element vector that establishes a
%                            nondefault spatial coordinate system by
%                            specifying the image XData. The value can
%                            have more than 2 elements, but only the first
%                            and last elements are actually used.
%
%   'YData'                  Two-element vector that establishes a
%                            nondefault spatial coordinate system by
%                            specifying the image YData. The value can
%                            have more than 2 elements, but only the first
%                            and last elements are actually used.
%
%   Class Support
%   -------------  
%   A truecolor image can be uint8, uint16, double, or single. An indexed image
%   can be logical, uint8, double, or single.  An intensity image can be
%   logical, uint8, double, single, int16, or uint16. The input image must be
%   nonsparse. 
%
%   If your image is int16 or single, the CData in the resulting image
%   object will be double. For all other classes, the CData matches the
%   input image class.
% 
%   Related Toolbox Preferences
%   ---------------------------  
%   You can use the IPTSETPREF function to set several toolbox preferences that
%   modify the behavior of IMSHOW:
%
%   - 'ImshowBorder' controls whether IMSHOW displays the image with a border
%     around it.
%
%   - 'ImshowAxesVisible' controls whether IMSHOW displays the image with the
%     axes box and tick labels.
%
%   - 'ImshowInitialMagnification' controls the initial magnification for
%     image display, unless you override it in a particular call by
%     specifying IMSHOW(...,'InitialMagnification',INITIAL_MAG).
%
%   For more information about these preferences, see the reference entry for
%   IPTSETPREF.
%   
%   Remarks
%   -------
%   IMSHOW is the toolbox's fundamental image display function, optimizing 
%   figure, axes, and image object property settings for image display. IMTOOL
%   provides all the image display capabilities of IMSHOW but also provides 
%   access to several other tools for navigating and exploring images, such as
%   the Pixel Region tool, Image Information tool, and the Adjust Contrast 
%   tool. IMTOOL presents an integrated environment for displaying images and
%   performing some common image processing tasks.
%
%   Examples
%   --------
%       % Display an image from a file
%       imshow('board.tif') 
%
%       % Display an indexed image
%       [X,map] = imread('trees.tif');
%       imshow(X,map) 
%
%       % Display an intensity image 
%       I = imread('cameraman.tif');
%       imshow(I) 
%
%       % Display an intensity image, adjust the diplay range
%       h = imshow(I,[0 80]);
%
%   See also IMREAD, IMTOOL, IPTGETPREF, IPTSETPREF, SUBIMAGE, TRUESIZE,
%            WARP, IMAGE, IMAGESC.

%   Copyright 1993-2004 The MathWorks, Inc.  
%   $Revision: 1.1.8.3 $  $Date: 2004/12/18 07:36:11 $

  varargin_translated = preParseInputs(varargin{:});
  
  [cdata, cdatamapping, clim, map, xdata, ydata, ...
   initial_mag, style] = ...
      imageDisplayParseInputs(varargin_translated{:});
    
  if isempty(initial_mag)
    initial_mag = iptgetpref('ImshowInitialMagnification');
  else
    initial_mag = checkInitialMagnification(initial_mag,{'fit'},...
                                            mfilename,'INITIAL_MAG', ...
                                            []);
  end
  
  new_figure = isempty(get(0,'CurrentFigure')) | ...
      strcmp(get(get(0,'CurrentFigure'), 'NextPlot'), 'new');
  
  if (new_figure)
      fig_handle = figure('Visible', 'off');
      ax_handle = axes('Parent', fig_handle);
  else
      ax_handle = newplot;
      fig_handle = ancestor(ax_handle,'figure');
  end
  

  do_fit = strcmp(initial_mag,'fit');
  if isempty(style)
    style = get(fig_handle,'WindowStyle');
  end
  if ~do_fit && strcmp(style,'docked')
      wid = sprintf('Images:%s:magnificationMustBeFitForDockedFigure',mfilename);
      warning(wid,'%s%s','The initial magnification of the image is set to',...
              ' ''fit'' in a docked figure.');
      do_fit = true;
  end
  
  hh = basicImageDisplay(fig_handle,ax_handle,...
                         cdata,cdatamapping,clim,map,xdata,ydata);
  set(get(ax_handle,'Title'),'Visible','on');
  set(get(ax_handle,'XLabel'),'Visible','on');
  set(get(ax_handle,'YLabel'),'Visible','on');
  
  single_image = isSingleImageDefaultPos(fig_handle, ax_handle);
  
  
  if single_image && do_fit && isBorderTight 
      % Have the image fill the figure.
      set(ax_handle, 'Units', 'normalized', 'Position', [0 0 1 1]);
      axes_moved = true;

  elseif single_image && ~do_fit 
      initSize(hh,initial_mag/100,isBorderTight)
      axes_moved = true;          

  else
      axes_moved = false;
      
  end

  if axes_moved
      % The next line is so that a subsequent plot(1:10) goes back
      % to the default axes position. 
      set(fig_handle, 'NextPlot', 'replacechildren');
  end      
  
  if (nargout > 0)
    % Only return handle if caller requested it.
    h = hh;
  end
  
  if (new_figure)
      set(fig_handle, 'Visible', 'on');
  end

end % imshow

%------------------------------------------------------------------------------
function is_tight = isBorderTight
  
     borderPref = iptgetpref('ImshowBorder');
     if strcmp(borderPref, 'tight')
         is_tight = true;
     else
         is_tight = false;
     end
     
end

%----------------------------------------------------------------------
function varargin_translated = preParseInputs(varargin)
% Catch old style syntaxes and warn

% Obsolete syntaxes:
%   IMSHOW(I,N) 
%   N is ignored, gray(256) is always used for viewing intensity images.
%
%   IMSHOW(...,DISPLAY_OPTION) 
%   DISPLAY_OPTION is translated as follows:
%   'truesize'   -> 'InitialMagnification', 100
%   'notruesize' -> 'InitialMagnification', 'fit'
%
%   IMSHOW(x,y,A,...) 
%   x,y are translated to 'XData',x,'YData',y

new_args = {};
num_args = nargin;

if (num_args == 0)
    eid = sprintf('Images:%s:tooFewArgs',mfilename);
    error(eid,'%s\n%s','IMSHOW expected at least 1 input argument',...
          'but was called instead with 0 input arguments.')
end

if (num_args > 1) && ischar(varargin{end}) 
    % IMSHOW(...,DISPLAY_OPTION)

    str = varargin{end};
    strs = {'truesize', 'notruesize'};
    try
        % If trailing string is not 'truesize' or 'notruesize' jump to
        % catch block and pass trailing string argument to regular input
        % parsing so error will come from that parsing code.
        option = iptcheckstrs(str, strs, mfilename,'DISPLAY_OPTION', nargin);

        % Remove old argument
        varargin(end) = [];  
        num_args = num_args - 1;
    
        % Translate to InitialMagnification
        new_args{1} = 'InitialMagnification';
        if strncmp(option,'truesize',length(option))
            new_args{2} = 100;
            msg1 = 'IMSHOW(...,''truesize'') is an obsolete syntax. ';
            msg2 = 'Use IMSHOW(...,''InitialMagnification'',100) instead.';
        
        else
            new_args{2} = 'fit';
            msg1 = 'IMSHOW(...,''notruesize'') is an obsolete syntax. ';
            msg2 = 'Use IMSHOW(...,''InitialMagnification'',''fit'') instead.';
        
        end

        wid = sprintf('Images:%s:obsoleteSyntaxDISPLAY_OPTION',mfilename);        
        warning(wid,'%s\n%s',msg1,msg2)
    catch
        % Trailing string did not match 'truesize' or 'notruesize' let regular
        % parsing deal with it.
        
        % Reset lasterr since we are ignoring the error from iptcheckstrs if 
        % a valid syntax was used.
        lasterr('');
    end
end

if (num_args==3 || num_args==4) && ...
            isvector(varargin{1}) && isvector(varargin{2}) && ...
            isnumeric(varargin{1}) && isnumeric(varargin{2})             
    % IMSHOW(x,y,...)

    % Translate to IMSHOW(...,'XData',x,'YData',y)
    p = length(new_args);
    new_args{p+1} = 'XData';
    new_args{p+2} = varargin{1};
    new_args{p+3} = 'YData';
    new_args{p+4} = varargin{2};
    
    % Remove old arguments
    varargin(1:2) = [];
    num_args = num_args - 2;

    wid = sprintf('Images:%s:obsoleteSyntaxXY',mfilename);            
    msg1 = 'IMSHOW(x,y,...) is an obsolete syntax. ';
    msg2 = 'Use IMSHOW(...,''XData'',x,''YData'',y) instead.';    
    warning(wid,'%s%s',msg1,msg2)
end

if num_args == 2 && (numel(varargin{2}) == 1)
    % IMSHOW(I,N)

    wid = sprintf('Images:%s:obsoleteSyntaxN',mfilename);                
    msg1 = 'IMSHOW(I,N) is an obsolete syntax. Your intensity ';
    msg2 = 'image will be displayed using 256 shades of gray.';
    warning(wid,'%s%s',msg1,msg2)

    % Remove old argument
    varargin(2) = [];
end

varargin_translated = {varargin{:}, new_args{:}};

end
