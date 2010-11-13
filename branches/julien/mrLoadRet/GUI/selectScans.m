function scanList = selectScans(view,title,groupNum)
% scanList = selectScans(view,[title]);
%
%   Gather a list of scans available in Inplane/TSeries
%   and query the user for a sub-selection.
%
%   An alternate functionchooseScans uses numScans(view)
%   to determine the number of scans to choose from
%   Use selectScans if you will be analyzing the tSeries. 
%   Use chooseScans, if your code does not depend on the 
%   presence/absence of the tSeries files.
%
% Output:
%  scanList: list of selected scans.
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

%Check for zero:
if nScans == 0
  mrErrorDlg('No scans found!');
  return
end

for i = 1:nScans
  scanNames{i} = sprintf('%i:%s (%s)',i,viewGet(view,'description',i,groupNum),viewGet(view,'tSeriesFile',i,groupNum));
end

% Which scans to analyze?
iSel = buttondlg(title, scanNames);
scanList = find(iSel);

return;
