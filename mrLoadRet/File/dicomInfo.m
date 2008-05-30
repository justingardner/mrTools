% dicomInfo.m
%
%        $Id$
%      usage: dicomInfo(scanNum,groupNum)
%         by: justin gardner
%       date: 05/02/07
%    purpose: print out dicom information from scan
%       e.g.: dicomInfo(1);
function retval = dicomInfo(scanNum,groupNum)

% check arguments
if ~any(nargin == [2])
  help dicomInfo
  return
end

v = newView;
dicom = viewGet(v,'dicom',scanNum,groupNum);
originalFileName = viewGet(v,'originalFileName',scanNum,groupNum);
originalGroupName = viewGet(v,'originalGroupName',scanNum,groupNum);

if isempty(dicom)
  disp(sprintf('(dicomInfo) No dicom info for scan: %i in group: %s',scanNum,viewGet(v,'groupName',groupNum)));
  return
end

% get info out of dicom
paramsInfo{1} = {'dicomFileNum',1,'round=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',length(dicom))};
dicomFields = fieldnames(dicom{1}.ACQ);
for i = 1:length(dicom)
  % put in originalFile/Group
  paramsInfo{2}{1} = 'original';
  if length(originalFileName)>=i
    paramsInfo{2}{2}{i} = sprintf('%s: %s',originalGroupName{i},originalFileName{i});
  else
    paramsInfo{2}{2}{i} = 'unknown';
  end    
  paramsInfo{2}{3} = 'editable=0';
  paramsInfo{2}{4} = 'contingent=dicomFileNum';
  paramsInfo{2}{5} = 'type=string';
  % now get each field of the dicom header
  for j = 1:length(dicomFields)
    paramsInfo{j+2}{1} = dicomFields{j};
    paramsInfo{j+2}{2}{i} = dicom{i}.ACQ.(dicomFields{j});
    paramsInfo{j+2}{3} = 'editable=0';
    paramsInfo{j+2}{4} = 'contingent=dicomFileNum';
    paramsInfo{j+2}{5} = 'type=string';
  end
end

mrParamsDialog(paramsInfo,'Dicom parameters');

  
  