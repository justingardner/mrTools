% makeEmptyROI.m
%
%        $Id:$ 
%      usage: makeEmptyROI(v,<'scanNum'>,<'groupNum'>,<'name=...'>)
%             makeEmptyROI(v,<'scanNum=...'>,<'groupNum=...'>,<'name=...'>)
%             makeEmptyROI(v,<'scanNum',scanNum>,<'groupNum',groupNum>,<'name',name>))
%         by: justin gardner
%       date: 12/31/11
%    purpose: creates an empty roi with coordinates set for the scan and group
%             if scanNum is not defined, the current scan of the current group is used
%             if scanNum = 0, the current base is used instead
%
function roi = makeEmptyROI(v,varargin)

% check arguments
if nargin < 1
  help makeEmptyROI
  return
end

if ~isview(v),disp(sprintf('(makeEmptyROI) Invalid view passed in'));return,end

% get arguments
scanNum = [];groupNum = [];
getArgs(varargin,{'scanNum=[]','groupNum=[]','name=[]'});

% make a name
if ieNotDefined('name')
  % go through roi names and get the largest numbered
  % roi name, i.e. ROI4 then make then new name ROI5
  maxnum = 0;
  roiNames = viewGet(v,'roiNames');
  for i = 1:length(roiNames)
    if regexp(roiNames{i},'^ROI\d+$')
      maxnum = max(maxnum,str2num(roiNames{i}(4:end)));
    end
  end
  name=sprintf('ROI%.0f',maxnum+1);
end
roi.name = name;
if isequal(scanNum,0)
  baseNum = viewGet(v,'currentBase');
  roi.sformCode = viewGet(v,'baseSformCode',baseNum);
  % if the baseSformCode == 0 then use the bases Qform matrix
  % since it means that the alignment has not been set.
  if viewGet(v,'baseSformCode',baseNum) == 0
    roi.xform = viewGet(v,'baseQform',baseNum);
    roi.sformCode = 0;
  else
    roi.xform = viewGet(v,'baseSform',baseNum);
  end
  roi.baseNum = viewGet(v,'currentBase');
  roi.voxelSize = viewGet(v,'baseVoxelSize',baseNum);
else
  roi.voxelSize = viewGet(v,'scanVoxelSize',scanNum,groupNum);
  if viewGet(v,'scanSformCode',scanNum,groupNum)
    roi.xform = viewGet(v,'scanSform',scanNum,groupNum);
  else
    roi.xform = viewGet(v,'scanQform',scanNum,groupNum);
  end
end

[tf roi] = isroi(roi);





