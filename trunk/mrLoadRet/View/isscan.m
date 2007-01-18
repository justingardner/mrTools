function val = isscan(scanParams)
% function val = isview(scanParams)
% 
% djh, 2007

val = 1;

if ieNotDefined('scanParams')
    val = 0;
    return
end

if ~isstruct(scanParams)
	val = 0;
	return
end

requiredFields = {'description','fileName','fileType','niftiHdr',...
		  'voxelSize','totalFrames','junkFrames','nFrames',...
		  'dataSize','framePeriod','originalFileName','originalGroupName'};

for f = 1:length(requiredFields)
  fieldName = requiredFields{f};
  if ~isfield(scanParams,fieldName)
    mrWarnDlg(['Invalid scanParams, missing field:',fieldName]);
    val = 0;
  end
end

if isfield(scanParams,'originalFileName') && isfield(scanParams,'originalGroupName')
  if length(scanParams.originalFileName) ~= length(scanParams.originalGroupName)
    mrWarnDlg('Invalid scanParams: originalFileName does not match originalGroupName');
  end
end
