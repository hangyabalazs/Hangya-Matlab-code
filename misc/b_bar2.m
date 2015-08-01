function [xo,yo] = bar2(color,varargin)
%BAR2 Colored bar graph.
%   The only difference between BAR and BAR2 is that using BAR2 you have to
%   specify color in the first input argument. See BAR for detailed description
%   of calling syntax.
%
%   See also BAR, BARH, HIST and PLOT.

error(nargchk(1,4,nargin));

[msg,x,y,xx,yy,linetype,plottype,barwidth,equal] = makebars(varargin{:});
if ~isempty(msg), error(msg); end

if nargout==2,
  warning(sprintf(...
     ['BAR with two output arguments is obsolete.  Use H = BAR(...) \n',...
      '         and get the XData and YData properties instead.']))
  xo = xx; yo = yy; % Do not plot; return result in xo and yo
else % Draw the bars
  cax = newplot;
  next = lower(get(cax,'NextPlot'));
  hold_state = ishold;
  edgec = get(gcf,'defaultaxesxcolor');
  facec = 'flat';
  h = []; 
  cc = ones(size(xx,1),1);
  if ~isempty(linetype), facec = linetype; end
  for i=1:size(xx,2)
    numBars = (size(xx,1)-1)/5;
    f = 1:(numBars*5);
    f(1:5:(numBars*5)) = [];
    f = reshape(f, 4, numBars);
    f = f';

    v = [xx(:,i) yy(:,i)];

    h=[h patch('faces', f, 'vertices', v, 'cdata', i*cc, ...
        'FaceColor',color,'EdgeColor',edgec)];
  end
  if length(h)==1, set(cax,'clim',[1 2]), end
  if ~equal, 
    hold on,
    plot(x(:,1),zeros(size(x,1),1),'*')
  end
  if ~hold_state, 
    % Set ticks if less than 16 integers
    if all(all(floor(x)==x)) & (size(x,1)<16),  
      set(cax,'xtick',x(:,1))
    end
    hold off, view(2), set(cax,'NextPlot',next);
    set(cax,'Layer','Bottom','box','on')
    % Turn off edges when they start to overwhelm the colors
    if size(xx,2)*numBars > 150, 
       set(h,{'edgecolor'},get(h,{'facecolor'}));
    end
  end
  if nargout==1, xo = h; end
end
