% getSmoothColor.m
%
%        $Id: getSmoothColor.m,v 1.2 2008/07/16 21:38:18 justin Exp $ 
%      usage: getSmoothColor(colorNum,totalColors,<colorMap>,<skipRange>)
%         by: justin gardner
%       date: 07/08/08
%    purpose: returns a color that smoothly varies over totalColors
%             colorNum is a number from 1 to totalColors specifying the desired color
%             totalColors is the total number of colors that can be returned
%             colorMap is the colormap function to use (e.g. 'hot', 'cool', 'hsv', ...)
%             if skipRange is specified, then you can prevent colors from going all
%             white, by specifying that you want to use only the bottom 80% of colors
%             for instance (i.e. by specifying skipRange to be 0.8); If you want
%             to use the top 80% of colors, make skipRange = -0.8;
%       e.g.: 
%
%for i = 1:10
%  plot(i,1,'ko','MarkerFaceColor',getSmoothColor(i,10),'MarkerSize',16);hold on
%  text(i,1,sprintf('%i',i),'HorizontalAlignment','Center','Color',getSmoothColor(11-i,10));
%end
%
function color = getSmoothColor(colorNum,totalColors,colorMap,skipRange)

% default return color
color = [0.5 0.5 0.5];

% check arguments
if ~any(nargin == [1 2 3 4])
  help getSmoothColor
  return
end

% default arguments
if ieNotDefined('totalColors'), totalColors = 256;end
if ieNotDefined('colorMap'), colorMap = 'gray';end
if ~any(strcmp(colorMap,{'hsv','gray','pink','cool','bone','copper','flag'}))
  if ~exist(colorMap,'file')
    disp(sprintf('(getSmoothColor) Unknown colormap function %s',colorMap));
    return
  end
end
if ieNotDefined('skipRange')
  if strcmp(colorMap,'gray')
    skipRange = 0.8;
  else
    skipRange = 1;
  end
end

% get colors to choose from
if skipRange > 0
  colors = eval(sprintf('%s(ceil(totalColors*((1-skipRange)+1)))',colorMap));
else
  colors = eval(sprintf('%s(ceil(totalColors*((1+skipRange)+1)))',colorMap));
  colors = colors(end-totalColors+1:end,:);
end  

% select out the right color
if (colorNum >= 1) & (colorNum <= totalColors)
  color = colors(colorNum,:);
else
  % out of bounds. Warn and return gray
  disp(sprintf('(getSmoothColor) Color %i out of bounds [1 %i]',colorNum,totalColors));
end


