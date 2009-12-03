% importCaret.m
%
%        $Id:$ 
%      usage: mlrImportCaret()
%         by: justin gardner
%       date: 12/03/09
%    purpose: function to import caret surfaces
%
function retval = mlrImportCaret(varargin)

% check arguments
if ~any(nargin == [0 1])
  help importCaret
  return
end

% we need files from these places
surfRelaxDir = 'surfRelax';
caretFileDir = 'surf';
atlasDir = '../PALS_B12.LR';

% check to make sure the directories exist
if ~checkDirs(caretFileDir,surfRelaxDir),return,end

% look for coord files in caret directory. These will be used for computing
% the xform that goes from 711-2B back to the original coordinates
[rightCoordFiles d.rightStemName] = getCoordFiles(caretFileDir,'R.Midthickness');
if isempty(rightCoordFiles),return,end
[leftCoordFiles d.leftStemName] = getCoordFiles(caretFileDir,'L.Midthickness');
if isempty(leftCoordFiles),return,end

% now compute the transformation from caret 711-2B back to original coordinates;
d.leftXform = computeTransformBetweenSurfaces(leftCoordFiles{2},leftCoordFiles{1});
if isempty(d.leftXform),return,end
d.rightXform = computeTransformBetweenSurfaces(rightCoordFiles{2},rightCoordFiles{1});
if isempty(d.rightXform),return,end

% now look for the files in the atlas directory
[topo d.leftAtlasAligned] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('matchStr=%s',d.leftStemName),'noDisplay=1');
if isempty(d.leftAtlasAligned),disp(sprintf('(mlrImportCaret) Could not find deformed surface for %s in %s',d.leftStemName,atlasDir));return,end
[topo d.rightAtlasAligned] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('matchStr=%s',d.rightStemName),'noDisplay=1');
if isempty(d.rightAtlasAligned),disp(sprintf('(mlrImportCaret) Could not find deformed surface for %s in %s',d.rightStemName,atlasDir));return,end

% look for topo files for LEFT and RIGHT hemispheres
[d.leftAtlasTopo d.rightAtlasTopo] = getAtlasTopo(atlasDir);
if isempty(d.leftAtlasTopo) || isempty(d.rightAtlasTopo),return,end

% get all the display surfaces that match
[topo d.leftDisplaySurface] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('numNodes=%i',d.leftAtlasTopo.numNodes),sprintf('matchStr=LEFT'),'noDisplay=1');
if isempty(d.leftDisplaySurface),disp(sprintf('(mlrImportCaret) Could not find any display coordinages for LEFT in atlas dir: %s',atlasDir)),return,end
[topo d.rightDisplaySurface] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('numNodes=%i',d.rightAtlasTopo.numNodes),sprintf('matchStr=RIGHT'),'noDisplay=1');
if isempty(d.rightDisplaySurface),disp(sprintf('(mlrImportCaret) Could not find any display coordinages for RIGHT in atlas dir: %s',atlasDir)),return,end

% get the volume header
[hdr d.volumeFileName] = getNiftiVolumeHeader(surfRelaxDir);

% get how much to shift coords to be in the center
d.coordShift = hdr.dim(2:4)/2;
  
% display the transforms
dispXform('Left hemisphere 711-2B to original xform',d.leftXform);
dispXform('Right hemisphere 711-2B to original xform',d.rightXform);

% Now load left surfaces
d.leftCoordSurf = loadSurfCaret(fullfile(atlasDir,d.leftAtlasAligned.name),fullfile(atlasDir,d.leftAtlasTopo.name),'xform',d.leftXform,'zeroBased=1','coordShift',d.coordShift);
d.leftFiducialSurf = loadSurfCaret(fullfile(atlasDir,d.leftDisplaySurface(1).name),fullfile(atlasDir,d.leftAtlasTopo.name));
d.leftInflatedSurf = loadSurfCaret(fullfile(atlasDir,d.leftDisplaySurface(2).name),fullfile(atlasDir,d.leftAtlasTopo.name));
d.leftVeryInflatedSurf = loadSurfCaret(fullfile(atlasDir,d.leftDisplaySurface(3).name),fullfile(atlasDir,d.leftAtlasTopo.name));

% Now load right surfaces
d.rightCoordSurf = loadSurfCaret(fullfile(atlasDir,d.rightAtlasAligned.name),fullfile(atlasDir,d.rightAtlasTopo.name),'xform',d.rightXform,'zeroBased=1','coordShift',d.coordShift);
d.rightFiducialSurf = loadSurfCaret(fullfile(atlasDir,d.rightDisplaySurface(1).name),fullfile(atlasDir,d.rightAtlasTopo.name));
d.rightInflatedSurf = loadSurfCaret(fullfile(atlasDir,d.rightDisplaySurface(2).name),fullfile(atlasDir,d.rightAtlasTopo.name));
d.rightVeryInflatedSurf = loadSurfCaret(fullfile(atlasDir,d.rightDisplaySurface(3).name),fullfile(atlasDir,d.rightAtlasTopo.name));

% compute curvature
doCurvature = 1;
if doCurvature
  d.leftAtlasCurvature = calcCurvature(d.leftFiducialSurf);
  d.rightAtlasCurvature = calcCurvature(d.rightFiducialSurf);
end

% check for save directory
if ~isdir('caret'),mkdir('caret');end

% and write out left surfaces
writeOFF(d.leftCoordSurf,fullfile('caret','leftCoords'));
writeOFF(d.leftFiducialSurf,fullfile('caret','leftAtlasFiducial'));
writeOFF(d.leftInflatedSurf,fullfile('caret','leftAtlasInflated'));
writeOFF(d.leftVeryInflatedSurf,fullfile('caret','leftAtlasVeryInflated'));
if doCurvature,saveVFF(fullfile('caret','leftAtlasCurvature'),d.leftAtlasCurvature);end

% write out right surfaces
writeOFF(d.rightCoordSurf,fullfile('caret','rightCoords'));
writeOFF(d.rightFiducialSurf,fullfile('caret','rightAtlasFiducial'));
writeOFF(d.rightInflatedSurf,fullfile('caret','rightAtlasInflated'));
writeOFF(d.rightVeryInflatedSurf,fullfile('caret','rightAtlasVeryInflated'));
if doCurvature,saveVFF(fullfile('caret','rightAtlasCurvature'),d.rightAtlasCurvature);end

% and link anatomy file
linkFile(d.volumeFileName,fullfile('caret',getLastDir(d.volumeFileName)));
linkFile(setext(d.volumeFileName,'img'),fullfile('caret',setext(getLastDir(d.volumeFileName),'img')));

%%%%%%%%%%%%%%%%%%%%%%
%%   getAtlasTopo   %%
%%%%%%%%%%%%%%%%%%%%%%
function [leftAtlasTopo rightAtlasTopo] = getAtlasTopo(atlasDir)

leftAtlasTopo = listCaretDir(sprintf('dirname=%s',atlasDir),'matchStr=LEFT','noDisplay=1');
if isempty(leftAtlasTopo)
  disp(sprintf('(mlrImportCaret) Could not find any atlas topo files that match LEFT in %s',atlasDir));
  return
end
rightAtlasTopo = listCaretDir(sprintf('dirname=%s',atlasDir),'matchStr=RIGHT','noDisplay=1');
if isempty(rightAtlasTopo)
  disp(sprintf('(mlrImportCaret) Could not find any atlas topo files that match RIGHT in %s',atlasDir));
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
  disp(sprintf('(mlrImportCaret) Could not find %s.coord file',matchStr));
  return
end
if length(coordFile) ~= 1
  disp(sprintf('(mlrImportCaret) Found %i files that match %s.coord file, but should have found only 1',matchStr));
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
  disp(sprintf('(mlrImportCaret) Could not find file %s',caretName));
  return
end

filenames{1} = fullfile(caretFileDir,coordFile.name);
filenames{2} = caretName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkAllCaretFiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkAllCaretFiles(surfRelaxDir,caretFileDir)

% check directory
disppercent(-inf,sprintf('(mlrImportCaret) Checking %s directory for topo and coord files',caretFileDir));
[topoFiles coordFiles] = listCaretDir(sprintf('dirname=%s',caretFileDir),'noDisplay=1');
disppercent(inf);

% return if files not found
if isempty(topoFiles)
  disp(sprintf('(mlrImportCaret) Could not find any topo files in %s. Have you run mypreborder/mypostborder?',caretFileDir));
  return
end
% return if files not found
if isempty(coordFiles)
  disp(sprintf('(mlrImportCaret) Could not find any coord files in %s. Have you run mypreborder/mypostborder?',caretFileDir));
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
  disp(sprintf('(mlrImportCaret) Could not find any nifti volume files in %s',surfRelaxDir));
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
  disp(sprintf('(mlrImportCaret) Num nodes between coord (%i) and topo (%i) do not match',coordFiles(params.coordFileNum).numNodes,topoFiles(params.topoFileNum).numNodes));
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
function [hdr volumeFileName] = getNiftiVolumeHeader(surfRelaxDir)

hdr = [];
volumeFileName = [];
% get nifti volume files
niftiVolumes = dir(fullfile(surfRelaxDir,'*.hdr'));
niftiVolumes = cat(1,niftiVolumes,dir(fullfile(surfRelaxDir,'*.nii')));
% check that we got something
if isempty(niftiVolumes)
  disp(sprintf('(mlrImportCaret:getNiftiVolumeHeader) Could not find any nifti volume files in %s',surfRelaxDir));
  return
end
if length(niftiVolumes) > 1
  paramsInfo{1} = {'niftiVolume',dir2cell(niftiVolumes),'The nifti volume should be the canonical volume from which this surface was made'};

  % bring up dialog box
  params = mrParamsDialog(paramsInfo,'Select topo and coord file to display');
  if isempty(params),return,end
  volumeFileName = fullfile(surfRelaxDir,params.niftiVolume);
  hdr = cbiReadNiftiHeader(volumeFileName);
else
  volumeFileName = fullfile(surfRelaxDir,niftiVolumes.name);
  hdr = cbiReadNiftiHeader(volumeFileName);
end


%%%%%%%%%%%%%%%%%%%
%%   checkDirs   %%
%%%%%%%%%%%%%%%%%%%
function isok = checkDirs(varargin)

isok = 0;
for i = 1:nargin
  if ~isdir(varargin{i})
    disp(sprintf('(mlrImportCaret) Could not find %s directory',varargin{i}));
    return
  end
end
isok = 1;

