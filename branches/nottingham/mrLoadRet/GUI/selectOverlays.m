function overlayList = selectOverlays(thisView,title)
% overlayList = selectOverlays(thisView,[title]);
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

overlayNames = viewGet(thisView,'overlayNames');

%Check for zero:
if isempty(overlayNames)
  mrWarnDlg('(selectOverlays) No overlay found!');
  return
end

% Which overlays to analyze?
iSel = buttondlg(title, overlayNames);
overlayList = find(iSel);

return;
