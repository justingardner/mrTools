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
iSel = buttondlg(title, baseNames,preselection);
if isempty(iSel)
  basesList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  basesList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
end

return;
