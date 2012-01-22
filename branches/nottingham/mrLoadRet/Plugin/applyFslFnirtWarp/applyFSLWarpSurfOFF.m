% applyFSLWarpSurfOFF.m.m
%
%       $Id: applyFSLWarpSurfOFF.m.m 1833 2010-11-13 18:37:14Z julien $	
%      usage: applyFSLWarpSurfOFF.m()
%         by: julien besle
%       date: 10/05/2011
%    purpose: applies FNIRT non-linear transformation to OFF surface files
%
%             assumes that the FNIRT input volume has an sform matrix that puts it in the canonical base space


function applyFSLWarpSurfOFF

startPathStr = pwd;
filterspec = {'*.img','FNIRT warp coeff .img file'; '*.*','All files'};
title = 'Choose structural to functional FNIRT coefficient file';
warpCoefFileName = getPathStrDialog(startPathStr,title,filterspec,'off');
if isempty(warpCoefFileName),return,end;

startPathStr = fileparts(warpCoefFileName);
filterspec = {'*.img','NIFTI .img file'; '*.*','All files'};
title = 'Choose FNIRT input volume';
inputVolumeFilename = getPathStrDialog(startPathStr,title,filterspec,'off');
if isempty(inputVolumeFilename),return,end

startPathStr = mrGetPref('volumeDirectory');
filterspec = {'*.img','NIFTI .img file'; '*.*','All files'};
title = 'Choose canonical base volume';
canonicalBaseFilename = getPathStrDialog(startPathStr,title,filterspec,'off');
if isempty(canonicalBaseFilename),return,end

startPathStr = fileparts(canonicalBaseFilename);
filterspec = {'*.off','SurfRelax OFF file';'*.off','SurfRelax off surface file'; '*.*','All files'};
title = 'Choose surface OFF file(s) to convert';
surfFileNames = getPathStrDialog(startPathStr,title,filterspec,'on');
% Aborted
if isempty(surfFileNames),return,end


canonicalHdr = cbiReadNiftiHeader(canonicalBaseFilename);
inputHdr = cbiReadNiftiHeader(inputVolumeFilename);

% load the surface coordinates
coordsCount = 0;
allCoords = [];
cSurf = 0;
isFlat = false(1,length(surfFileNames));
for iSurf = 1:length(surfFileNames)
  cSurf = cSurf+1;
  surf{iSurf} = loadSurfOFF(surfFileNames{iSurf});
  if isfield(surf{iSurf},'patch2parent')
    isFlat(iSurf)=true;
  else
    allCoords = [allCoords;surf{iSurf}.vtcs];
    coordsIndices(iSurf,:) = coordsCount+[1 surf{iSurf}.Nvtcs];
    coordsCount = coordsCount+surf{iSurf}.Nvtcs;
  end
end
    

if coordsCount
  %make coordinate matrix homogenous
  allCoords = [allCoords ones(coordsCount,1)]';

  %convert coordinates into indices in the input volume

  % use jonas' function to get the array2world xform
  array2world = mlrXFormFromHeader(canonicalHdr,'array2world');

  %convert coordinates from world to canonical base array space
  allCoords = array2world\allCoords;
  allCoords = shiftOriginXform*allCoords;
  %from base array space to canonical base space
  allCoords = canonicalHdr.sform44*allCoords;
  %from  canonical base space to input volume array space
  allCoords = inputHdr.sform44\allCoords;
  allCoords = shiftOriginXform\allCoords;


  warpedCoords = applyFSLWarpCoords(allCoords,inputHdr.pixdim(2:4)',[.5 .5 .5], warpCoefFileName, 'temp.img', inputHdr, 1);

  warpedCoords = shiftOriginXform*warpedCoords;
  warpedCoords = inputHdr.sform44*warpedCoords;
  warpedCoords = canonicalHdr.sform44\warpedCoords;
  warpedCoords = shiftOriginXform\warpedCoords;
  warpedCoords = array2world*warpedCoords;
end


%first save surfaces
for iSurf = find(~isFlat)
  surf{iSurf}.vtcs = warpedCoords(1:3,coordsIndices(iSurf,1):coordsIndices(iSurf,2))';
  [pathname,filename,extension] = fileparts(surfFileNames{iSurf});
  outputFilename =  [pathname,'/',filename,'_invFNIRT',extension];
  [filename,pathname]=uiputfile('*.off',['Save surface converted from ' filename],outputFilename);
  newFilename{iSurf} = fullfile(pathname,filename);
  if ~isnumeric(filename)
    writeOFF(surf{iSurf}, newFilename{iSurf});
  end
end

%now save flat files with new parent file name (it's actually the only thing that changes)
%in fact it's kind of useless because when importing the file in mrLoad it is possible to change the parent surface...
%but at least I learnt something
for iSurf = find(isFlat)
  saveFile = true;
  [isconverted,whichSurface] = ismember(surf{iSurf}.parentSurfaceName,surfFileNames);
  [pathname,filename,extension] = fileparts(surfFileNames{iSurf});
  if ~isconverted
    surf{iSurf}.parentSurfaceName = getPathStrDialog(startPathStr,['Choose new parent surface (WM) for flat file converted from ' filename],{'*WM*.off','Inner surface OFF file'});
    if isempty(surf{iSurf}.parentSurfaceName)
      saveFile=false;
    end
  else
    saveFile=true;
    surf{iSurf}.parentSurfaceName = newFilename{whichSurface};
  end
  if saveFile
    outputFilename =  [pathname,'/',filename,'_invFNIRT',extension];
    [filename,pathname]=uiputfile('*.off',['Save flat surface converted from ' filename],outputFilename);
    outputFilename = fullfile(pathname,filename);
    if ~isnumeric(filename)
      writePatchOFF(surf{iSurf}, outputFilename)
    end
  end
end  


% writeOFF.m
%
%      usage: writeOFF()
%         by: eli merriam, taken and adapted from makeFlat by julien besle
%       date: 11/05/2011
%    purpose: 
%
function retval = writePatchOFF(surf, fileName)

% check arguments
if ~any(nargin == [2])
  help writeOFF
  return
end

% write the OFF format file 
fid = fopen(fileName, 'w', 'ieee-be');
fprintf(fid, '#PATCH\n');
fprintf(fid, '#parent_surface=%s\n', surf.parentSurfaceName);
fprintf(fid, '#parent_dimensions=%i %i %i\n', surf.nParent(1), surf.nParent(2), surf.nParent(3));
fprintf(fid, '#patch_dimensions=%i %i %i\n', surf.Nvtcs, surf.Ntris, 1 );
fprintf(fid, '#parent_vertex_indexes:\n');
for i=1:length(surf.patch2parent)
  fprintf(fid, '#%i %i\n', surf.patch2parent(i,1)-1, surf.patch2parent(i,2)-1);
end

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [size(surf.vtcs,1) size(surf.tris,1) 0], 'int32'); 

% Vertices
fwrite(fid, surf.vtcs', 'float32');

% Faces
fwrite(fid, [3*ones(size(surf.tris,1),1) surf.tris-1 zeros(size(surf.tris,1),1)]', 'int32');

% Close file
fclose(fid);

return;

  
  
  
