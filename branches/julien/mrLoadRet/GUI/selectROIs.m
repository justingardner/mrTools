function roisList = selectROIs(thisView,title,preselected)
% roisList = selectROIs(thisView,[title]);
%
%   Gather a list of loaded ROIs 
%   and query the user for a sub-selection.
%
% Output:
%  roisList: list of selected ROIs.
%
% 06/10/10 jb adapted from selectScans
%
% $Id$ 

roisList = [];

if ieNotDefined('title')
  title = 'Choose ROIs';
end
if ieNotDefined('preselected')
  preselected = [];
end

roiNames = viewGet(thisView,'roiNames');

%Check for zero:
if isempty(roiNames)
  mrWarnDlg('(selectROIs) No ROI found!');
  return
end

preselection = zeros(1,length(roiNames));
preselection(preselected) = 1;
iSel = buttondlg(title, roiNames,preselection);

if isempty(iSel)
  roisList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  roisList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
end

return;
