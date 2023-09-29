% mlrGetSurfaceNames.m
%
%      usage: [leftSurfaces rightSurfaces] = mlrGetSurfaceNames()
%         by: justin gardner
%       date: 09/05/19
%    purpose: uses mrSurfViewer to have user select surfaces
%
function [leftSurfaces, rightSurfaces] = mlrGetSurfaceNames(varargin)

getArgs(varargin,{'surfPath=[]'});

% now find surfaces if we weren't passed in the names
if isempty(surfPath)
  titleStr = 'Choose left gray matter (outer) surface for this subject';
  disp(sprintf('(mlrImportFreesurferLAbel) %s',titleStr));
  surfPath = mlrGetPathStrDialog(mrGetPref('volumeDirectory'),titleStr,'*.off');
  if isempty(surfPath)
    leftSurfaces = [];
    rightSurfaces = [];
    return;
  end
end

% make sure this is a GM (not WM or Inf)
surfFilename = getLastDir(surfPath);
surfPath = fileparts(surfPath);
surfFilename = replaceStr(surfFilename,'WM.','GM.');
surfFilename = replaceStr(surfFilename,'Inf.','GM.');

% Make names for left surfaces
leftSurfFilename = replaceStr(surfFilename,'right','left');
leftSurfaces.outerSurface = replaceStr(leftSurfFilename,'GM.','Inf.');
leftSurfaces.outerCoords = leftSurfFilename;
leftSurfaces.innerSurface = replaceStr(leftSurfFilename,'GM.','WM.');
leftSurfaces.innerCoords = replaceStr(leftSurfFilename,'GM.','WM.');
leftSurfaces.curv = setext(replaceStr(leftSurfFilename,'GM.','Curv.'),'vff');
leftSurfaces.path = surfPath;

% Make names for right surfaces
rightSurfFilename = replaceStr(surfFilename,'left','right');
rightSurfaces.outerSurface = replaceStr(rightSurfFilename,'GM.','Inf.');
rightSurfaces.outerCoords = rightSurfFilename;
rightSurfaces.innerSurface = replaceStr(rightSurfFilename,'GM.','WM.');
rightSurfaces.innerCoords = replaceStr(rightSurfFilename,'GM.','WM.');
rightSurfaces.curv = setext(replaceStr(rightSurfFilename,'GM.','Curv.'),'vff');
rightSurfaces.path = surfPath;

% now look for antomy file
surfDir = dir(surfPath);
anatNames = {};
for iFile = 1:length(surfDir)
  if any(strcmp(getext(surfDir(iFile).name),{'hdr','nii'}))
    anatNames{end+1} = surfDir(iFile).name;
  end
end

% no anatomy files
if length(anatNames) < 1
  disp(sprintf('(mlrGetSurfaceNames) Could not find 3D anatomy in folder: %s',surfPath));
  return
end

% too many anatomy files
if length(anatNames) > 1
  paramsInfo{1} = {'anatomyName',anatNames};
  params = mrParamsDialog(paramsInfo,'Choose anatomy file');
  if isempty(params),return,end
  anatNames = {params.anatomyName};
end

% set surface names
leftSurfaces.anatomy = anatNames{1};
rightSurfaces.anatomy = anatNames{1};


%%%%%%%%%%%%%%%%%%%%
%    replaceStr    %
%%%%%%%%%%%%%%%%%%%%
function s = replaceStr(s,searchstr,replacestr)

foundloc = strfind(s,searchstr);
if ~isempty(foundloc)
  s = sprintf('%s%s%s',s(1:(foundloc(1)-1)),replacestr,s((foundloc(1)+length(searchstr)):end));
end

