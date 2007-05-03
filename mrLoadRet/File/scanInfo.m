% scanInfo.m
%
%      usage: scanInfo(scanNum,groupNum,<displayInDialog>)
%         by: justin gardner
%       date: 11/08/06
%    purpose: print out information about scans in group
%       e.g.: groupInfo(1);
%             groupInfo('Raw');
function retval = scanInfo(scanNum,groupNum,displayInDialog)

% check arguments
if ~any(nargin == [2 3])
  help scanInfo
  return
end

if ~exist('displayInDialog','var'),displayInDialog = 0;end

view = newView('Volume');

% if groupNum is a string, then user passed in a name rather than
% a number
if isstr(groupNum)
  groupName = groupNum;
  groupNum = viewGet(view,'groupNum',groupName);
  if isempty(groupNum)
    disp(sprintf('(groupInfo): No group %s',groupName));
    return
  end
end

if (groupNum < 1) | (groupNum > viewGet(view,'numberOfGroups'))
  disp(sprintf('(groupInfo): No group number %i',groupNum));
  return
end
groupName = viewGet(view,'groupName',groupNum);
% set the group and scan
view = viewSet(view,'curGroup',groupNum);

% grab info
description = viewGet(view,'description',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);
tr = viewGet(view,'framePeriod',scanNum,groupNum);
totalFrames = viewGet(view,'totalFrames',scanNum,groupNum);
filename = viewGet(view,'tSeriesFile',scanNum,groupNum);
originalFilename = viewGet(view,'originalFilename',scanNum,groupNum);
originalGroupname = viewGet(view,'originalGroupname',scanNum,groupNum);
stimFilename = viewGet(view,'stimFilename',scanNum,groupNum);
scanHdr = viewGet(view,'niftiHdr',scanNum,groupNum);
scanDims = viewGet(view,'scanDims',scanNum,groupNum);

% get params info
matfile = viewGet(view,'params',scanNum,groupNum);
fieldsToIgnore = {'tseriesfiles','groupname','description','filename'};
% display info
if displayInDialog
  paramsInfo = {{'description',description,'editable=0','Scan description'},...
		{'Filename',filename,'editable=0','Name of file'},...
		{'GroupName',groupName,'editable=0','Name of group'}};
  for i = 1:length(originalFilename)
    paramsInfo{end+1} = {sprintf('Original%i',i) sprintf('%s: %s',originalGroupname{i},originalFilename{i}) 'editable=0' 'Name of original group and filename that this scan came from'};
  end
  for i = 1:length(stimFilename)
    paramsInfo{end+1} = {sprintf('Stimfile%i',i) stimFilename{i} 'editable=0' 'Name of stim file'};
  end
  paramsInfo{end+1} = {'voxelSize',scanVoxelSize,'editable=0','Voxel dimensions in mm'};
  paramsInfo{end+1} =  {'dims',scanDims,'editable=0','Dimensions of scan'};
  paramsInfo{end+1} =  {'TR',tr,'editable=0','TR'};
  paramsInfo{end+1} =  {'numVolumes',totalFrames,'editable=0','Number of volumes'};
  paramsInfo{end+1} = {'qform',scanHdr.qform44,'editable=0','Qform matrix specifies the transformation to the scanner coordinate frame'};
  paramsInfo{end+1} = {'sform',scanHdr.sform44,'editable=0','Sform matrix is set by mrAlign and usually specifies the transformation to the volume anatomy'};
  paramsInfo{end+1} = {'sformCode',scanHdr.sform_code,'editable=0','If sformCode is 0 it means the sform has never been set and mrLoadRet will use the qform to compute the transform to the base anatomy. If mrAlign has been run properly, then this value should be set to 1'};

  % display parameters form associated matfile
  if ~isempty(matfile) && isfield(matfile,'params')
    fields = fieldnames(matfile.params);
    for i = 1:length(fields)
      % ignore fields from the list set above
      if ~any(strcmp(lower(fields{i}),fieldsToIgnore))
	% and fields that have groupname in them
	if isempty(strfind(lower(fields{i}),'groupname'))
	  % display fields that are numeric
	  if isnumeric(matfile.params.(fields{i}))
	    paramsInfo{end+1} = {fields{i} num2str(matfile.params.(fields{i})) 'editable=0','Parameter from associated mat file'};
	  % or are strings
	  elseif isstr(matfile.params.(fields{i}))
	    paramsInfo{end+1} = {fields{i} matfile.params.(fields{i}) 'editable=0','Parameter from associated mat file'};
	  end
	end
      end
    end
  end
  
  mrParamsDialog(paramsInfo,'Scan info');
else
  disp(sprintf('%s',description));
  disp(sprintf('Filename: %s GroupName: %s',filename,groupName));
  for i = 1:length(originalFilename)
    disp(sprintf('Original Filename: %s Group: %s',originalFilename{i},originalGroupname{i}));
  end
  for i = 1:length(stimFilename)
    disp(sprintf('StimFilename: %s',stimFilename{i}));
  end
  
  disp(sprintf('voxelSize=[%0.1f %0.1f %0.1f] TR=%0.4f Dims: [%i %i %i] Volumes=%i',scanVoxelSize(1),scanVoxelSize(2),scanVoxelSize(3),tr,scanDims(1),scanDims(2),scanDims(3),totalFrames));

  % display qform and sform
  disp(sprintf('++++++++++++++++++++++++++ qform ++++++++++++++++++++++++++'));
  for rownum = 1:4
    disp(sprintf('%f\t%f\t%f\t%f',scanHdr.qform44(rownum,1),scanHdr.qform44(rownum,2),scanHdr.qform44(rownum,3),scanHdr.qform44(rownum,4)));
  end

  disp(sprintf('++++++++++++++++++++++++++ sform ++++++++++++++++++++++++++'));
  for rownum = 1:4
    disp(sprintf('%f\t%f\t%f\t%f',scanHdr.sform44(rownum,1),scanHdr.sform44(rownum,2),scanHdr.sform44(rownum,3),scanHdr.sform44(rownum,4)));
  end
end


