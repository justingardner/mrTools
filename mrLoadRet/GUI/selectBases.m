function basesList = selectBases(thisView,title,preselected)
% basesList = selectROIs(thisView,[title]);
%
%   Gather a list of loaded Bases 
%   and query the user for a sub-selection.
%
% Output:
%  baseList: list of selected ROIs.
%
% 28/10/10 jb adapted from selectScans
%
% $Id$ 

basesList = [];

if ieNotDefined('title')
  title = 'Choose Bases';
end
if ieNotDefined('preselected')
  preselected = [];
end

baseNames = viewGet(thisView,'baseNames');

%Check for zero:
if isempty(baseNames)
  mrWarnDlg('(selectROIs) No ROI found!');
  return
end

preselection = zeros(1,length(baseNames));
preselection(preselected) = 1;
basesList = find(buttondlg(title, baseNames,preselection));

return;
