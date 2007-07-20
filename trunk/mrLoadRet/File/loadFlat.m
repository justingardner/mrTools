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

% Strip extension to make sure it is .mat
pathStr = [stripext(pathStr),'.mat'];

% File does not exist
if ~exist(pathStr,'file')
  mrWarnDlg(['File ',pathStr,' not found']);
  return
end

% load the flat file
flat = load(pathStr);

% check its fields
if any(~isfield(flat,{'curvature','gLocs2d','gLocs3d'}))
  mrWarnDlg(sprintf('(loadFlat) %s is not a flat file generated using SurfRelax',pathStr));
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now load the base anatomy file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filterspec = {'*.hdr','Nifti file header'};
title = 'Choose volume anatomy file header for this flat image';
anatPathStr = getPathStrDialog(fileparts(pathStr),title,filterspec);

% Strip extension to make sure it is .mat
anatPathStr = [stripext(anatPathStr),'.hdr'];

% File does not exist
if ~exist(anatPathStr,'file')
  mrWarnDlg(['File ',anatPathStr,' not found']);
  return
end

% load the anatomy file header
hdr = cbiReadNiftiHeader(anatPathStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the flat image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first get coordinates
xmin = min(flat.gLocs2d(:,1));
xmax = max(flat.gLocs2d(:,1));
ymin = min(flat.gLocs2d(:,2));
ymax = max(flat.gLocs2d(:,2));
x = xmin:xmax;
y = ymin:ymax;

% now interp the curvature
disppercent(-inf,'Creating flat image');
for i = 1:length(x)
  disppercent(i/length(x));
  for j = 1:length(y)
    % find nearest point in curvature
    dist = (flat.gLocs2d(:,1)-x(i)).^2+(flat.gLocs2d(:,2)-y(j)).^2;
    flat.pos(i,j) = first(find(min(dist)==dist));
    flat.map(i,j,1) = flat.curvature(flat.pos(i,j));
    flat.map(i,j,2) = flat.curvature(flat.pos(i,j));
    flat.baseCoord(i,j,:) = flat.gLocs3d(flat.pos(i,j),:);
  end
end
disppercent(inf);

% Extract permutation matrix to keep track of slice orientation.
% This logic which is admittedly arcane is duplicated in mrAlignGUI. If you
% make changes here, please update that function as well.
[q,r] = qr(inv(hdr.qform44(1:3,1:3)));
permutationMatrix = abs([q(1,:); q(2,:); q(3,:)]);


% now generate a base 
base.data = flat.map;
base.hdr = hdr;
base.name = getLastDir(pathStr);
base.permutationMatrix = permutationMatrix;
base.coords = flat.baseCoord;


v = viewSet(v,'newBase',base);


