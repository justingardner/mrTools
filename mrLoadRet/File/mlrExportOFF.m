% mlrExportOff.m
%
%        $Id:$ 
%      usage: mlrExportOff(filename,off,<vertexColors>)
%         by: justin gardner
%       date: 12/23/18
%    purpose: Function to convert off surface files into output formats
%             compatible with 3D printing and importing into rendering
%             programs - STL and OFF
%
%             filename is the file to save under - use extension obj for
%             waverfront obj format or stl for stl format
%
%             off is a structure with the fields vtcs and tris or
%             an off filename
%
%             vertexColors is an optional matrix that is nx1 or nx3
%             where n is the number of vertices and can be either
%             grayscale color of vertex (0-1) or color [R G B] also
%             in values of 0 - 1.
% 
%             e.g.
%             mlrExportOff('ouptut.obj','jg_left.off');
%
%             vertexColors can also be a vff file:
%             mlrExportOff('jg_left_GM.obj','jg_left_GM.off','jg_left_Curv.vff');
%
%
function retval = mlrExportOFF(filename,off,vertexColors)

% check arguments
retval = [];
if ~any(nargin == [2 3])
  help mlrExportOFF
  return
end
if nargin < 3
  vertexColors = [];
end

% check to see if off is a string in which case it should be loaded as a file
if isstr(off)
  if ~mlrIsFile(off)
    disp(sprintf('(mlrExportOFF) OFF surface file %s not found',off));
    return
  else
    % load the surface
    off = loadSurfOFF(off);
    % if empty then return, bad file
    if isempty(off),return,end
  end
end

% make sure off is a off sturcture
if ~isstruct(off) || ~all(isfield(off,{'vtcs','tris'}))
  disp(sprintf('(mlrExportOFF) Passed in structure is not a off returned by loadSurfOFF'));
  return
end

% get the filetype from the extension
filetype = lower(getext(filename));
if isempty(filetype)
  % default to obj
  filetype = 'obj';
  filename = setext(filename,filetype);
end

% get name of surface from filename
[~,surfaceName] = fileparts(filename);
surfaceName = stripext(surfaceName);

% make sure we have a recognized extension
if ~any(strcmp(filetype,{'obj','stl'}))
  disp(sprintf('(mlrExportOFF) Unrecognized file format %s. Should be stl or obj',filetype));
  return
end

% now open file for saving
outfile = fopen(filename,'w');
if outfile == -1
  disp(sprintf('(mlrExportOFF) Could not open output file %s',filename));
  return
end

% check vertexColors
if ~isempty(vertexColors)
  % check if it is a filename
  if isstr(vertexColors)
    if isfile(vertexColors)
      [data hdr] = loadVFF(vertexColors);
      if isempty(data),return,end
      vertexColors = data(:);
    else
      disp(sprintf('(mlrExportOFF) %s is not a vff file',vertexColors));
      return
    end
  end
  if strcmp(filetype,'stl')
    disp(sprintf('(mlrExportOFF) STL format cannot handle vertex colors. Ignoring'));
    vertexColors = [];
  end
  if size(vertexColors,1) ~= size(off.vtcs,1)
    disp(sprintf('(mlrExportOFF) Vertex colors must have one entry for every vertex. Ignoring\n'));
    vertexColors = [];
  end
  if ~any(size(vertexColors,3) == [1 3])
    disp(sprintf('(mlrExportOFF) Vertex colors must be a grayscale color a r,g,b triplet. Ignoring\n'));
    vertexColors =  [];
  end
end

if strcmp(filetype,'stl')
  writeSTL(outfile,off.vtcs,off.tris,surfaceName);
elseif strcmp(filetype,'obj')
  writeOBJ(outfile,off.vtcs,off.tris,vertexColors);
end

% close the file
fclose(outfile);

%%%%%%%%%%%%%%%%%%
%    writeOBJ    %
%%%%%%%%%%%%%%%%%%
function writeOBJ(outfile,vtcs,tris,vertexColors)

% compute sizes
nVtcs = size(vtcs,1);
nTris = size(tris,1);

% reverse x-axis (otherwise the surfaces seem to flip)
vtcs(:,1) = 256-vtcs(:,1);

% write vertex colors
if ~isempty(vertexColors)
  if size(vertexColors,2) == 1
    % write vertices with grayscale
    disppercent(-inf,'(mlrExportOFF) Writing vertices with grayscale');
    for iVertex = 1:nVtcs
      fprintf(outfile,'v %f %f %f %f\n',vtcs(iVertex,1),vtcs(iVertex,2),vtcs(iVertex,3),vertexColors(iVertex,1));
      if mod(iVertex,5000) == 1,disppercent(iVertex/nVtcs);end
    end
    disppercent(inf);
  else
    % write vertices with rgb
    disppercent(-inf,'(mlrExportOFF) Writing vertices with color');
    for iVertex = 1:nVtcs
      fprintf(outfile,'v %f %f %f %f %f %f\n',vtcs(iVertex,1),vtcs(iVertex,2),vtcs(iVertex,3),vertexColors(iVertex,1),vertexColors(iVertex,2),vertexColors(iVertex,3));
      if mod(iVertex,5000) == 1,disppercent(iVertex/nVtcs);end
    end
    disppercent(inf);
  end    
else
  % write vertices with no colors
  disppercent(-inf,'(mlrExportOFF) Writing vertices');
  for iVertex = 1:nVtcs
    fprintf(outfile,'v %f %f %f\n',vtcs(iVertex,1),vtcs(iVertex,2),vtcs(iVertex,3));
    if mod(iVertex,5000) == 1,disppercent(iVertex/nVtcs);end
  end
  disppercent(inf);
end

% write faces
disppercent(-inf,'(mlrExportOFF) Writing faces');
for iTri = 1:nTris
  fprintf(outfile,'f %i %i %i\n',tris(iTri,1),tris(iTri,2),tris(iTri,3));
  if mod(iTri,5000) == 1,disppercent(iTri/nTris);end
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%
%    writeSTL    %
%%%%%%%%%%%%%%%%%%
function writeSTL(outfile,vtcs,tris,surfaceName)

% write out start of file
fprintf(outfile,'solid %s\n',surfaceName);


disppercent(-inf,'(mlrExportOFF) Converting triangles');
nTris = size(tris,1);
for iTri = 1:nTris
  % get vertices
  v1 = vtcs(tris(iTri,1),:);
  v2 = vtcs(tris(iTri,2),:);
  v3 = vtcs(tris(iTri,3),:);
  % compute normal
  n = cross(v2-v1,v3-v2);
  n = n/norm(n);
  % now print it to file
  fprintf(outfile,'facet normal %f %f %f\n',n(1),n(2),n(3));
  fprintf(outfile,'    outer loop\n');
  fprintf(outfile,'        vertex %f %f %f\n',v1(1),v1(2),v1(3));
  fprintf(outfile,'        vertex %f %f %f\n',v2(1),v2(2),v2(3));
  fprintf(outfile,'        vertex %f %f %f\n',v3(1),v3(2),v3(3));
  fprintf(outfile,'    endloop\n');
  fprintf(outfile,'endfacet\n');
  disppercent(iTri/nTris);
end
disppercent(inf);

fprintf(outfile,'endsolid %s\n',surfaceName);


