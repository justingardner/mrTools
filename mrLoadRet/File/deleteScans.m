% deleteScans.m
%
%        $Id$
%      usage: deleteScans(<v>,<scanList>,<groupName/Num>)
%         by: justin gardner
%       date: 09/19/07
%    purpose: delete scans
%
function v = deleteScans(v,scanList,group)

% check arguments
if ~any(nargin == [1 2 3])
  help mrDeleteScans
  return
end

% get a view if not passed in one.
if ieNotDefined(v)
  v = newView;
end

% set group if passed in
if ~ieNotDefined('group')
  v = viewSet(v,'curGroup',group);
end

% get a scan list if it has not been passed in
if ieNotDefined('scanList')
  scanList = selectScans(v);
end

% check the scan list
badScanNums = setdiff(scanList,1:viewGet(v,'nScans'));
if ~isempty(badScanNums)
  disp(sprintf('(mrDeleteScans) No scans %s',num2str(badScanNums)));
  return
end

for iScan = 1:length(scanList)
    % first get time series name for each one of these scans
    % since as we delete them, then numbers stop making sense
    tSeriesFile{iScan} = viewGet(v,'tSeriesFile',scanList(iScan));
end
% now go through and delete
for iScan = 1:length(scanList)
    % get the scan number
    scanNum = viewGet(v,'scanNum',tSeriesFile{iScan});
    if ~isempty(scanNum)
        scanNum = scanNum(1);
        v = viewSet(v,'deleteScan',scanNum);
        disp(sprintf('Scan for file %s deleted.',tSeriesFile{iScan}));
    else
        disp(sprintf('(mrDeleteScans) Could not delete scan for file %s',tSeriesFile{iScan}));
    end
    if ~isempty(v.figure)
      refreshMLRDisplay(v.viewNum);
    end
end
if ~isempty(scanList)
    disp(sprintf('To remove the nifti files for these deleted scans run mrCleanDir'));
end
