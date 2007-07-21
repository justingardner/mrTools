% loadFlat.m
%
%        $Id$
%      usage: loadFlat(v,flatFileName)
%         by: justin gardner
%       date: 07/20/07
%    purpose: 
%
function v = loadFlat(v,flatFileName)

% check arguments
if ~any(nargin == [0 1 2])
  help loadFlat
  return
end

% get mrGlobals and view
mrGlobals;

% Open dialog box to have user choose the file
if ieNotDefined('flatFileName')
  if ieNotDefined('flatFilePath')
    startPathStr = mrGetPref('volumeDirectory');
  else
    startPathStr = flatFilePath;
  end
  filterspec = {'*.mat','Matlab flat file'};
  title = 'Choose flat file';
  pathStr = getPathStrDialog(startPathStr,title,filterspec);
else
  pathStr = flatFileName;
end

% Aborted
if ieNotDefined('pathStr')
  return
end

% make into a cell array
pathStr = cellArray(pathStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now load the base anatomy file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterspec = {'*.hdr','Nifti file header'};
title = 'Choose volume anatomy file header for this flat image';
anatPathStr = getPathStrDialog(fileparts(pathStr{1}),title,filterspec);

% Strip extension to make sure it is .mat
anatPathStr = [stripext(anatPathStr),'.hdr'];

% File does not exist
if ~exist(anatPathStr,'file')
  mrWarnDlg(['File ',anatPathStr,' not found']);
  return
end

% load the anatomy file header
hdr = cbiReadNiftiHeader(anatPathStr);

% now go through all flat file names, load and convert
% into images
flat = {};maxWidth = -inf;maxHeight = -inf;
for fileNum = 1:length(pathStr)
  % Strip extension to make sure it is .mat
  filename = [stripext(pathStr{fileNum}),'.mat'];

  % File does not exist
  if ~exist(filename,'file')
    mrWarnDlg(['File ',filename,' not found']);
    continue
  end

  % load the flat file
  flat{end+1} = load(filename);

  % check its fields
  if ~isfield(flat{end},'curvature') || ~isfield(flat{end},'gLocs2d') || ~isfield(flat{end},'gLocs3d')
    mrWarnDlg(sprintf('(loadFlat{End}) %s is not a flat{end} file generated using SurfRelax',pathStr));
    continue
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % generate the flat{end} image
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % first get coordinates
  xmin = min(flat{end}.gLocs2d(:,1));
  xmax = max(flat{end}.gLocs2d(:,1));
  ymin = min(flat{end}.gLocs2d(:,2));
  ymax = max(flat{end}.gLocs2d(:,2));
  x = xmin:xmax;
  y = ymin:ymax;

  % now interp the curvature
  disppercent(-inf,sprintf('Creating flat image for %s',filename));
  for i = 1:length(x)
    disppercent(i/length(x));
    for j = 1:length(y)
      % find nearest point in curvature
      dist = (flat{end}.gLocs2d(:,1)-x(i)).^2+(flat{end}.gLocs2d(:,2)-y(j)).^2;
      flat{end}.pos(i,j) = first(find(min(dist)==dist));
      % points that are greater than a distance of 5 away are
      % probably not in the flat{end} patch so mask them out
      if (min(dist) < 5)
	flat{end}.mask(i,j) = 1;
	flat{end}.baseCoords(i,j,:) = flat{end}.gLocs3d(flat{end}.pos(i,j),:);
	flat{end}.map(i,j) = flat{end}.curvature(flat{end}.pos(i,j));
      else
	flat{end}.mask(i,j) = 0;
	flat{end}.baseCoords(i,j,:) = [0 0 0];
	flat{end}.map(i,j) = 0;
      end
    end
  end
  disppercent(inf);

  % now blur/upsample image
  flat{end}.blurMap(:,:) = blur(flat{end}.map(:,:));
  flat{end}.median = median(flat{end}.blurMap(:));
  flat{end}.thresholdMap(:,:) = (flat{end}.blurMap(:,:)>median(flat{end}.blurMap(:)))*0.5+0.25;
  %flat{end}.thresholdMap = blur(flat{end}.thresholdMap);
  flat{end}.thresholdMap(~flat{end}.mask(:)) = 0;
  flat{end}.blurMap(~flat{end}.mask(:)) = 0;

  % set the maximum width and height of the images
  maxWidth = max(size(flat{end}.thresholdMap,1),maxWidth);
  maxHeight = max(size(flat{end}.thresholdMap,1),maxHeight);
end

% no flats made, then return
if isempty(flat)
  return
end

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(hdr.qform44(1:3,1:3)));
permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);

% now generate a base structure
base.hdr = hdr;
base.name = getLastDir(pathStr{1});
base.clip = [0 1];
base.permutationMatrix = permutationMatrix;

% load all the flat maps into the base. We
% need to make all the flat images have
% the same width and height.
for i = 1:length(flat)
  base.data(:,:,i) = flat{i}.thresholdMap;
  base.coordMap.coords(:,:,i,1) = flat{end}.baseCoords(:,:,3);
  base.coordMap.coords(:,:,i,2) = hdr.dim(3)-flat{end}.baseCoords(:,:,1);
  base.coordMap.coords(:,:,i,3) = hdr.dim(4)-flat{end}.baseCoords(:,:,2);
  base.coordMap.dims = hdr.dim(2:4)';
end

v = viewSet(v,'newBase',base);


