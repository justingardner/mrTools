function overlayList = selectOverlays(thisView,title,preselected)
% overlayList = selectOverlays(thisView,[title],[preselected]);
%
%   Gather a list of overlays available in current group/analysis
%   and query the user for a sub-selection.
%
% Output:
%  overlayList: list of selected overlays.
%
% 05/10/10 jb adapted from selectScans
%
% $Id$ 

overlayList = [];

if ieNotDefined('title')
  title = 'Choose overlays';
end
if ieNotDefined('preselected')
  preselected = [];
end


overlayNames = viewGet(thisView,'overlayNames');

%Check for zero:
if isempty(overlayNames)
  mrWarnDlg('(selectOverlays) No overlay found!');
  return
end

preselection = zeros(1,length(overlayNames));
preselection(preselected) = 1;

% Which overlays to analyze?
iSel = buttondlg(title, overlayNames,preselection);
if isempty(iSel)
  overlayList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  overlayList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
end

return;
