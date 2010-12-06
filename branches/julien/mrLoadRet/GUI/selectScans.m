function [scanList,scanNames] = selectScans(view,title,groupNum,preselected)
% scanList = selectScans(view,[title],[groupNum],[preselected]);
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

if ieNotDefined('title')
  title = 'Choose scans';
end

if ieNotDefined('groupNum')
   groupNum = viewGet(view,'currentGroup');
end
nScans = viewGet(view,'nScans',groupNum);

if ieNotDefined('preselected')
  preselected = [];
end

%Check for zero:
if nScans == 0
  mrErrorDlg('No scans found!');
  return
end

for i = 1:nScans
  scanNames{i} = sprintf('%i:%s (%s)',i,viewGet(view,'description',i,groupNum),viewGet(view,'tSeriesFile',i,groupNum));
end

preselection = zeros(1,length(scanNames));
preselection(preselected) = 1;

% Which scans to analyze?
iSel = buttondlg(title, scanNames,preselection);
scanList = find(iSel);

return;
