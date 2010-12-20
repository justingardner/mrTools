function scanList = selectScans(thisView,varargin)
% scanList = selectScans(thisView,[title],[preselected]);
%
%   this function is deprecated, use scanList = selectInList(thisView,'scans',title,preselected)
%
%   Gather a list of scans available in Inplane/TSeries
%   and query the user for a sub-selection.
%
%
% Output:
%  scanList: list of selected scans.
%  scanNames: names of all the scans the user had to choose from
%
% 4/16/99  dbr Initial code
% 3/30/2001, djh, added optional title string
% 11/9/06 jlg mrLoadRet 4 conversion
%
% $Id$	

scanList = selectInList(thisView,'scans',varargin);
mrWarnDlg('(selectScans) selectScans is deprecated. Please use ''selectInList(view,''scans'',...)'' instead.');
return

if ieNotDefined('title')
  title = 'Choose scans';
end

if ieNotDefined('groupNum')
   groupNum = viewGet(thisView,'currentGroup');
end
nScans = viewGet(thisView,'nScans',groupNum);

if ieNotDefined('preselected')
  preselected = [];
end

%Check for zero:
if nScans == 0
  mrErrorDlg('No scans found!');
  return
end

for i = 1:nScans
  scanNames{i} = sprintf('%i:%s (%s)',i,viewGet(thisView,'description',i,groupNum),viewGet(thisView,'tSeriesFile',i,groupNum));
end

preselection = zeros(1,length(scanNames));
preselection(preselected) = 1;

% Which scans to analyze?
iSel = buttondlg(title, scanNames,preselection);
if isempty(iSel)
  scanList = iSel; %if cancel has been pressed, this will be a 0*0 matrix, 
  %but if the top close button has been pressed, it will be a 0*1 matrix
else
  scanList = find(iSel); %if OK is pressed but nothing has been selected, this will output a 1*0 array
end

return;
