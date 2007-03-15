function reply = buttondlg(headerStr,optionStr)
% function reply = buttondlg(headerStr,optionStr)
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
%4/11/98 gmb   Wrote it
% 15/07/02  fwc     it will now make multiple columns of options if the number is large
% 2002.07.24  rfd & fwc: fixed minor bug where short checkbox labels were invisible.

OptionsPerColumn=30; % max number of options in one column

if nargin<2
    disp('Error: "bottondlg" requires two inputs');
    return
end

if isunix
    fontSize = 10;
else
    fontSize = 9;
end

if strcmp(class(optionStr),'cell')
    optionStr = char(optionStr);
end


nOptions = size(optionStr,1);

ncols=ceil(nOptions/OptionsPerColumn);
nOptionsPerColumn=ceil(nOptions/ncols);

%scale factors for x and y axis coordinates
xs = 1.8;  
ys = 1.8;

%default sizes
butWidth=10;
botMargin = 0.2;
%height = nOptions+2+botMargin*2+.5;
height = nOptionsPerColumn+2+botMargin*2+.5;
% If we don't have a minimum colwidth (5 in this case), short strings don't show up.
colwidth=max(size(optionStr,2),5);%
% width = max(size(optionStr,2),length(headerStr))+2;
width = max(ncols*colwidth,length(headerStr))+2;
width = max(width,2*butWidth+2);

%open the figure
h = figure('MenuBar','none',...
    'Units','char',...
    'Resize','off',...
    'NumberTitle','off');

% set the figure position
figpos = get(h,'Position');
figpos(3:4) = [width*xs,height*ys];
global MLR;
if isfield(MLR,'figpos') && isfield(MLR.figpos,'buttondlg')
  figpos(1:2) = MLR.figpos.buttondlg(1:2);
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
        'Units','char',...
        'String',optionStr(optionNum,:),...
        'BackgroundColor',bkColor,...
        'Position',[x*xs+c*colwidth*xs,y*ys,colwidth*xs+(c+1)*colwidth*xs,ys],...
        'HorizontalAlignment','left',...
        'FontSize',fontSize);
	if optionNum>=(c+1)*nOptionsPerColumn
        c=c+1;
        y=y0;
    end
end

%Display the OK/Cancel buttons
x=1;
y=botMargin;
uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','char',...
    'Position',[x*xs,y*ys,butWidth*xs,ys],...
    'CallBack','uiresume',...
    'FontSize',fontSize,...
    'UserData','Cancel');

x = width-butWidth-1;
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
else
    reply = [];
end

figpos = get(h,'Position');
MLR.figpos.buttondlg = figpos;
close(h)

