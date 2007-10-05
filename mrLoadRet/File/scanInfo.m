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

if (scanNum < 1) || (scanNum > viewGet(view,'nScans'))
  disp(sprintf('(scanInfo) Could not find scan %i in group %s',scanNum,groupName));
  return
end

% grab info
disppercent(-inf,'(scanInfo) Gathering scan info');
description = viewGet(view,'description',scanNum,groupNum);
scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);
tr = viewGet(view,'TR',scanNum,groupNum);
framePeriod = viewGet(view,'framePeriod',scanNum,groupNum);
totalFrames = viewGet(view,'totalFrames',scanNum,groupNum);
filename = viewGet(view,'tSeriesFile',scanNum,groupNum);
originalFilename = viewGet(view,'originalFilename',scanNum,groupNum);
originalGroupname = viewGet(view,'originalGroupname',scanNum,groupNum);
stimFilename = viewGet(view,'stimFilename',scanNum,groupNum);
scanHdr = viewGet(view,'niftiHdr',scanNum,groupNum);
scanDims = viewGet(view,'scanDims',scanNum,groupNum);
junkFrames = viewGet(view,'junkFrames',scanNum,groupNum);
totalJunkedFrames = viewGet(view,'totalJunkedFrames',scanNum,groupNum);

% get params info
nCols = 1;
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
  if ~isempty(stimFilename)
    % make stimfiles into sets of 4 (otherwise it takes too much room in the
    % dialog for long concatenations or averages).
    for j = 1:ceil(length(stimFilename)/4)
      stimfileNames = getLastDir(stimFilename{(j-1)*4+1});
      firstStimfile = ((j-1)*4+1);
      lastStimfile = min(length(stimFilename),j*4);
      for i = firstStimfile+1:lastStimfile
	stimfileNames = sprintf('%s, %s',stimfileNames,getLastDir(stimFilename{i}));
      end
      if lastStimfile == firstStimfile
	paramsInfo{end+1} = {sprintf('stimFilenames%i',lastStimfile),stimfileNames,'editable=0','Names of associated stimfiles'};
      else
	paramsInfo{end+1} = {sprintf('stimFilenames%ito%i',firstStimfile,lastStimfile),stimfileNames,'editable=0','Names of associated stimfiles'};
      end
    end
  end
  paramsInfo{end+1} = {'voxelSize',scanVoxelSize,'editable=0','Voxel dimensions in mm'};
  paramsInfo{end+1} =  {'dims',scanDims,'editable=0','Dimensions of scan'};
  paramsInfo{end+1} =  {'TR',tr,'editable=0','TR. This is retrieved frome the dicom header information.'};
  paramsInfo{end+1} =  {'framePeriod',framePeriod,'editable=0','Time each volume takes to acquire. Note that this is usually the same as TR, except for 3D scans in which this should be number of slices times the TR.'};
  paramsInfo{end+1} =  {'numVolumes',totalFrames,'editable=0','Number of volumes'};
  paramsInfo{end+1} =  {'junkFrames',junkFrames,'editable=0','Number of junk frames'};
  paramsInfo{end+1} =  {'totalJunkedFrames',num2str(totalJunkedFrames),'type=string','editable=0','Number of junk frames that have already been discarded from this time series (this is useful for aligning the volume number of a stimulus set in a stimulus file with the volumes in the time series)'};
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
  disppercent(inf);
  mrParamsDialog(paramsInfo,'Scan info',nCols);
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


