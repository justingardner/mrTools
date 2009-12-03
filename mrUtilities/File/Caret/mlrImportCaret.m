% importCaret.m
%
%        $Id:$ 
%      usage: importCaret()
%         by: justin gardner
%       date: 12/03/09
%    purpose: function to import caret surfaces
%
function retval = importCaret()

% check arguments
if ~any(nargin == [0 1])
  help importCaret
  return
end

% we need files from these places
surfRelaxDir = 'surfRelax';
caretFileDir = 'surf';

% check to make sure the directories exist
if ~isdir(caretFileDir)
  disp(sprintf('(importCaret) Could not find %s directory',caretFileDir));
  return
end

if ~isdir(surfRelaxDir)
  disp(sprintf('(importCaret) Could not find %s directory. Have you run mlrImportFreeSurfer?',surfRelaxDir));
  return
end

% look for right coord file
[rightCoordFiles leftStemName] = getCoordFiles(caretFileDir,'R.Midthickness');
if isempty(rightCoordFiles),return,end
[leftCoordFiles rightStemName] = getCoordFiles(caretFileDir,'L.Midthickness');
if isempty(leftCoordFiles),return,end

% display the names
disp(sprintf('(importCaret) Right hemisphere coords files are %s and %s',rightCoordFiles{1},rightCoordFiles{2}));
disp(sprintf('(importCaret) Left hemisphere coords files are %s and %s',leftCoordFiles{1},leftCoordFiles{2}));

% now compute the transformation from caret 711-2B back to original coordinates;
leftXform = computeTransformBetweenSurfaces(leftCoordFiles{2},leftCoordFiles{1});
if isempty(leftXform),return,end
rightXform = computeTransformBetweenSurfaces(rightCoordFiles{2},rightCoordFiles{1});
if isempty(rightXform),return,end

% display the transforms
dispXform('Left hemisphere 711-2B to original xform',leftXform);
dispXform('Right hemisphere 711-2B to original xform',rightXform);

% now look for the files in the atlas directory
atlasDir = '../PALS_B12.LR';
[topo leftAtlasAligned] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('matchStr=%s',leftStemName),'noDisplay=1');
if isempty(leftAtlasAligned),disp(sprintf('(importCaret) Could not find deformed surface for %s in %s',leftStemName,atlasDir));return,end
[topo rightAtlasAligned] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('matchStr=%s',rightStemName),'noDisplay=1');
if isempty(rightAtlasAligned),disp(sprintf('(importCaret) Could not find deformed surface for %s in %s',rightStemName,atlasDir));return,end

% look for topo files for LEFT and RIGHT hemispheres
leftAtlasTopo = listCaretDir(sprintf('dirname=%s',atlasDir),'matchStr=LEFT','noDisplay=1');
if isempty(leftAtlasTopo)
  disp(sprintf('(importCaret) Could not find any atlas topo files that match LEFT in %s',atlasDir));
  return
end
rightAtlasTopo = listCaretDir(sprintf('dirname=%s',atlasDir),'matchStr=RIGHT','noDisplay=1');
if isempty(rightAtlasTopo)
  disp(sprintf('(importCaret) Could not find any atlas topo files that match RIGHT in %s',atlasDir));
  return
end

% ask user to select topo files
paramsInfo{1} = {'leftTopo',dir2cell(leftAtlasTopo),'Choose the topo file for the left hemisphere'};
paramsInfo{2} = {'rightTopo',dir2cell(rightAtlasTopo),'Choose the topo file for the left hemisphere'};
params = mrParamsDialog(paramsInfo,'Select atlas topo files');
if isempty(params),return,end

% get the topo file info
for i = 1:length(leftAtlasTopo)
  if isequal(leftAtlasTopo(i).name,params.leftTopo)
    leftAtlasTopo = leftAtlasTopo(i);
    break;
  end
end
% get the topo file info
for i = 1:length(rightAtlasTopo)
  if isequal(rightAtlasTopo(i).name,params.rightTopo)
    rightAtlasTopo = rightAtlasTopo(i);
    break;
  end
end

% get all the display surfaces that match
[topo leftDisplaySurface] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('numNodes=%i',leftAtlasTopo.numNodes),sprintf('matchStr=LEFT'),'noDisplay=1');
if isempty(leftDisplaySurface),disp(sprintf('(importCaret) Could not find any display coordinages for LEFT in atlas dir: %s',atlasDir)),return,end
[topo rightDisplaySurface] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('numNodes=%i',rightAtlasTopo.numNodes),sprintf('matchStr=RIGHT'),'noDisplay=1');
if isempty(rightDisplaySurface),disp(sprintf('(importCaret) Could not find any display coordinages for RIGHT in atlas dir: %s',atlasDir)),return,end

% get the volume header
hdr = getNiftiVolumeHeader(surfRelaxDir);

% get how much to shift coords to be in the center
coordShift = hdr.dim(2:4)/2;
  
% Now load surfaces
% fix comments in loadSurCaret. Check for coordinate transform being correct
surf = loadSurfCaret(fullfile(atlasDir,leftAtlasAligned.name),fullfile(atlasDir,leftAtlasTopo.name),'xform',leftXform,'zeroBased=1','coordShift',coordShift);

writeOFF(surf,'test.off');
m = rand(1,surf.Nvtcs);
saveVFF('test.vff',m);

%%%%%%%%%%%%%%%%%%%
%    dispXform    %
%%%%%%%%%%%%%%%%%%%
function dispXform(name,xform)

disp(sprintf('(importCaret:dispXform) %s',name));
for i = 1:size(xform,1)
  disp(sprintf('%s',mynum2str(xform(i,:),'sigfigs=4')));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeTransformBetweenSurfaces    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xform = computeTransformBetweenSurfaces(fromSurfaceFilename,toSurfaceFilename)

xform = [];
% load the surfaces
fromSurface = openCaretFile(fromSurfaceFilename);
if isempty(fromSurface),return,end
toSurface = openCaretFile(toSurfaceFilename);
if isempty(toSurface),return,end

% check matching nodes
if fromSurface.num_nodes ~= toSurface.num_nodes
  disp(sprintf('(importCaret:computeTransformBetweenSurfaces) Nodes must match between 711-2B surface (%s:%i) and original surface (%s:%i)',fromSurfaceFilename,fromSurface.num_nodes,toSurfaceFilename,toSurface.num_nodes));
  return
end

% now compute the best matching xform between these two surfaces
toVertices = toSurface.data';toVertices(4,:) = 1;
fromVertices = fromSurface.data';fromVertices(4,:) = 1;
xform = toVertices*pinv(fromVertices);
xform(4,:) = [0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%
%    getCoordFiles    %
%%%%%%%%%%%%%%%%%%%%%%%
function [filenames stemName] = getCoordFiles(caretFileDir,matchStr)

filenames = [];
[topoFile coordFile] = listCaretDir(sprintf('dirname=%s',caretFileDir),sprintf('matchStr=%s.coord',matchStr),'noDisplay=1');
if isempty(coordFile)
  disp(sprintf('(importCaret) Could not find %s.coord file',matchStr));
  return
end
if length(coordFile) ~= 1
  disp(sprintf('(importCaret) Found %i files that match %s.coord file, but should have found only 1',matchStr));
  return
end  

% now go through and construct rest of files and make sure they are there
matchLoc = findstr(coordFile.name,matchStr);
stemName = sprintf('%s%s',coordFile.name(1:matchLoc-1),matchStr);

% construct the name of the surface which has been transformed into
% the Caret 711-2B space.
caretName = sprintf('%s_711-2B.coord',stemName);
caretName = fullfile(caretFileDir,caretName);
if ~isfile(caretName)
  disp(sprintf('(importCaret) Could not find file %s',caretName));
  return
end

filenames{1} = fullfile(caretFileDir,coordFile.name);
filenames{2} = caretName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkAllCaretFiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkAllCaretFiles(surfRelaxDir,caretFileDir)

% check directory
disppercent(-inf,sprintf('(importCaret) Checking %s directory for topo and coord files',caretFileDir));
[topoFiles coordFiles] = listCaretDir(sprintf('dirname=%s',caretFileDir),'noDisplay=1');
disppercent(inf);

% return if files not found
if isempty(topoFiles)
  disp(sprintf('(importCaret) Could not find any topo files in %s. Have you run mypreborder/mypostborder?',caretFileDir));
  return
end
% return if files not found
if isempty(coordFiles)
  disp(sprintf('(importCaret) Could not find any coord files in %s. Have you run mypreborder/mypostborder?',caretFileDir));
  return
end

% get topoFileInfo
for i = 1:length(topoFiles)
  topoFileNodes{i} = topoFiles(i).numNodes;
  topoFileComment{i} = topoFiles(i).comment;
end

% get coordFileInfo
for i = 1:length(coordFiles)
  coordFileNodes{i} = coordFiles(i).numNodes;
  coordFileComment{i} = coordFiles(i).comment;
end

% get nifti volume files
niftiVolumes = dir(fullfile(surfRelaxDir,'*.hdr'));
niftiVolumes = cat(1,niftiVolumes,dir(fullfile(surfRelaxDir,'*.nii')));

% check that we got something
if isempty(niftiVolumes)
  disp(sprintf('(importCaret) Could not find any nifti volume files in %s',surfRelaxDir));
  return
end

paramsInfo{1} = {'topoFileNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',length(topoFiles)),'Topo files contain the topology, i.e. the list of triangles that tile the surface'};
paramsInfo{end+1} = {'topoFileName',dir2cell(topoFiles),'type=String','group=topoFileNum','editable=0','Topo filename'};
paramsInfo{end+1} = {'topoFileNodes',topoFileNodes,'group=topoFileNum','type=numeric','editable=0','Number of nodes in topo file. Must match coord file.'};
paramsInfo{end+1} = {'topoFileComment',topoFileComment,'group=topoFileNum','type=string','editable=0','Comment for the topo files'};
paramsInfo{end+1} = {'coordFileNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',length(coordFiles)),'Coord files contain the coordinates of each vertex in the surface'};
paramsInfo{end+1} = {'coordFileName',dir2cell(coordFiles),'type=String','group=coordFileNum','editable=0','Name of coord file'};
paramsInfo{end+1} = {'coordFileNodes',coordFileNodes,'type=numeric','group=coordFileNum','editable=0','Number of nodes in coord file. Must match topo file.'};
paramsInfo{end+1} = {'coordFileComment',coordFileComment,'type=String','group=coordFileNum','editable=0','Comment for the coord file'};
paramsInfo{end+1} = {'niftiVolume',dir2cell(niftiVolumes),'The nifti volume should be the canonical volume from which this surface was made'};

% bring up dialog box
params = mrParamsDialog(paramsInfo,'Select topo and coord file to display');
if isempty(params),return,end

% check matching nodes
if topoFiles(params.topoFileNum).numNodes ~= coordFiles(params.coordFileNum).numNodes
  disp(sprintf('(importCaret) Num nodes between coord (%i) and topo (%i) do not match',coordFiles(params.coordFileNum).numNodes,topoFiles(params.topoFileNum).numNodes));
  return
end

% get the files
topoFileName = fullfile(caretFileDir,params.topoFileName{params.topoFileNum});
coordFileName = fullfile(surfRelaxDir,params.coordFileName{params.coordFileNum});
niftiFileName = fullfile(surfRelaxDir,params.niftiVolume);

% load the surface
surf = loadSurfCaret(coordFileName,topoFileName);
keyboard

%%%%%%%%%%%%%%%%%%
%    dir2cell    %
%%%%%%%%%%%%%%%%%%
function c = dir2cell(d)
		    
for i = 1:length(d)
  c{i} = d(i).name;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getNiftiVolumeHeader    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = getNiftiVolumeHeader(surfRelaxDir)
hdr = [];
% get nifti volume files
niftiVolumes = dir(fullfile(surfRelaxDir,'*.hdr'));
niftiVolumes = cat(1,niftiVolumes,dir(fullfile(surfRelaxDir,'*.nii')));
% check that we got something
if isempty(niftiVolumes)
  disp(sprintf('(importCaret) Could not find any nifti volume files in %s',surfRelaxDir));
  return
end
if length(niftiVolumes) > 1
  paramsInfo{1} = {'niftiVolume',dir2cell(niftiVolumes),'The nifti volume should be the canonical volume from which this surface was made'};

  % bring up dialog box
  params = mrParamsDialog(paramsInfo,'Select topo and coord file to display');
  if isempty(params),return,end
  hdr = cbiReadNiftiHeader(fullfile(surfRelaxDir,params.niftiVolume));
else
  hdr = cbiReadNiftiHeader(fullfile(surfRelaxDir,niftiVolumes.name));
end


