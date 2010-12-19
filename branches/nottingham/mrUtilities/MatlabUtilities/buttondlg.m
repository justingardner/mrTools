function reply = buttondlg(headerStr,optionStr,defaultValue)
% function reply = buttondlg(headerStr,optionStr,<defaultValue>)
%
% Button Dialog Box
%
%  Allows the user to toggle and select options in
%  the cell array string 'optionStr'.  Reply is a
%  boolean vector with length(optionStr) that indicates
%  selected options. reply is empty if the user
%  chooses to cancel.
%
%  Example:
%  reply = buttondlg('pick it',{'this','that','the other'})
%
%  Example:
%  reply = buttondlg('pick it',{'this','that','the other'},[0 1 0])
%
%4/11/98 gmb   Wrote it
% 15/07/02  fwc     it will now make multiple columns of options if the number is large
% 2002.07.24  rfd & fwc: fixed minor bug where short checkbox labels were invisible.
% 2010/08/27 jb adapt widht of each column to max length string

OptionsPerColumn=25; % max number of options in one column

if nargin<2
    disp('Error: "bottondlg" requires two inputs');
    return
end
if nargin<3
  defaultValue = zeros(1,length(optionStr));
else
  if length(defaultValue) < length(optionStr)
    defaultValue(end+1:length(optionStr)) = 1;
  end
end
  
  

if isunix
    fontSize = 10;
else
    fontSize = 9;
end

if ischar(optionStr)
   optionStr = cellstr(optionStr);
end

nOptions = length(optionStr);

ncols=ceil(nOptions/OptionsPerColumn);
nOptionsPerColumn=ceil(nOptions/ncols);

%scale factors for x and y axis coordinates
xs = 1.2;
ys = 1.5;

%default sizes
butWidth=10;
botMargin = 0.2;
%height = nOptions+2+botMargin*2+.5;
height = nOptionsPerColumn+2+botMargin*2+.5;

%compute a width per column
for iColumn=1:ncols
   %convert to char to find max length more easily
   tempOptionString = char(optionStr((iColumn-1)*nOptionsPerColumn+1:min(iColumn*nOptionsPerColumn,nOptions)));
   % If we don't have a minimum colwidth (5 in this case), short strings don't show up.
   colwidth(iColumn)=max(size(tempOptionString,2),5)+2;%
end
% width = max(size(optionStr,2),length(headerStr))+2;
width = max(sum(colwidth),length(headerStr))+2;
width = max(width,2*butWidth+2);

%open the figure
h = figure('MenuBar','none',...
    'Units','char',...
    'Resize','off',...
    'NumberTitle','off');


% set the figure position
figpos = get(h,'Position');
figpos(3:4) = [width*xs,height*ys];

figloc = mrGetFigLoc('buttondlg');
if ~isempty(figloc)
    figpos(1:2) = figloc(1:2);
end
set(h,'Position',figpos);

bkColor = get(h,'Color');

%Display title text inside a frame
x = 1;
% y = nOptions+1+botMargin;
y = nOptionsPerColumn+1+botMargin;

uicontrol('Style','frame',...
    'Units','char',...
    'String',headerStr,...
    'Position',[x*xs,(y+.2)*ys,(width-2)*xs,ys*1.3],...
    'HorizontalAlignment','center',...
    'FontSize',fontSize);

uicontrol('Style','text',...
    'Units','char',...
    'String',headerStr,...
    'Position',[(x+.25)*xs,(y+.3)*ys,(width-2.5)*xs,ys*.9],...
    'HorizontalAlignment','center',...
    'FontSize',fontSize);

%Display the radio buttons
y = y+botMargin/2;
y0=y;
c=0;

for optionNum=1:nOptions
    y = y-1;
    h_button(optionNum) = ...
        uicontrol('Style','checkbox',...
        'Value',defaultValue(optionNum),...
        'Units','char',...
        'String',optionStr{optionNum},...
        'BackgroundColor',bkColor,...
        'Position',[x*xs+sum(colwidth(1:c))*xs+c,y*ys,colwidth(c+1)*xs,ys],...
        'HorizontalAlignment','left',...
        'FontSize',fontSize);
    if optionNum>=(c+1)*nOptionsPerColumn
        c=c+1;
        y=y0;
    end
end

%Display the OK/Cancel buttons
%x=1;
x=min(.4*width-butWidth/2,.5*width-butWidth)*xs;
y=botMargin;
uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','char',...
    'Position',[x*xs,y*ys,butWidth*xs,ys],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','Cancel');

%x = width-butWidth-1;
x=max(.6*width-butWidth/2,.5*width)*xs;
uicontrol('Style','pushbutton',...
    'String','OK',...
    'Units','char',...
    'Position',[x*xs,y*ys,butWidth*xs,ys],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','OK');

%let the user select some radio buttons and
%wait for a 'uiresume' callback from OK/Cancel
uiwait

%determine which button was hit.
response = get(gco,'UserData');

%gather the status of the radio buttons if 'OK' was
%selected.  Otherwise return empty matrix.
if strcmp(response,'OK')
  for optionNum=1:nOptions
      reply(optionNum)=get(h_button(optionNum),'Value');
  end
elseif strcmp(response,'Cancel') %here we differentiate between a Cancel press
    reply = [];
else                             %and the user closing the window with the left (right) top button
    reply = zeros(0,1);
end

if ishandle(h)
  figpos = get(h,'Position');
  mrSetFigLoc('buttondlg',figpos);
  close(h)
  drawnow;
end
