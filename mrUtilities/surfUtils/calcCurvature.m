% calcCurvature.m
%
%        $Id$ 
%      usage: calcCurvature(innerSurfaceName,outerSurfaceName,<'nSmooth=2'>,<'vertexList',1:nVtcs>)
%         by: eli merriam & justin gardner
%       date: 08/05/08
%    purpose: compute curvature of the innerSurface. Pass in the names of the inner
%             and outer surfaces and it will return an array that has the curvature
%             values. You can save this using saveVFF. nSmooth is how many smoothing
%             passes to do. Default is 2. vertexList is the vertices you want
%             to compute the curvature on, default is all vertices.
%
function m = calcCurvature(innerSurfaceName,outerSurfaceName,varargin)

% check arguments
if ~any(nargin == [2 3])
  help calcCurvature
  return
end

nSmooth = [];vertexList = [];
getArgs(varargin,{'nSmooth=2','vertexList=[]'});

% off surface
innerSurfaceName = sprintf('%s.off',stripext(innerSurfaceName));
outerSurfaceName = sprintf('%s.off',stripext(outerSurfaceName));

% load the surfaces
if isfile(innerSurfaceName)
  innerSurf = loadSurfOFF(innerSurfaceName);
else
  disp(sprintf('(calcCurvature) Could not find file %s',innerSurfaceName));
  return
end
if isfile(outerSurfaceName)
  outerSurf = loadSurfOFF(outerSurfaceName);
else
  disp(sprintf('(calcCurvature) Could not find file %s',innerSurfaceName));
  return
end

% get the connection matrix between each vertex
surf.uniqueVertices = innerSurf.vtcs;
surf.uniqueFaceIndexList = innerSurf.tris;
connectionMatrix = findConnectionMatrix(surf);
clear surf;

connectionMatrix = double(connectionMatrix);
% compute vertex normals
vertexNormals = outerSurf.vtcs-innerSurf.vtcs;

% allocate space for m
m = zeros(1,innerSurf.Nvtcs);

disppercent(-inf,'(calcCurvature) Calculating curvature');
if isempty(vertexList),vertexList = 1:innerSurf.Nvtcs;end
for iVertex = vertexList
  % find neighbors of this vertex.
  vertexNeighbors = find(connectionMatrix(:,vertexList(iVertex)));
  % and also the neighbors of those neighbors
  [vertexNeighbors whichColumn] = find(connectionMatrix(:,vertexNeighbors));
  vertexNeighbors = unique(vertexNeighbors);
  [vertexNeighbors whichColumn] = find(connectionMatrix(:,vertexNeighbors));
  vertexNeighbors = unique(vertexNeighbors);
  % get x, y and z coordinates of neighbors relative to current vertex
  x = innerSurf.vtcs(vertexNeighbors,1)-innerSurf.vtcs(vertexList(iVertex),1);
  y = innerSurf.vtcs(vertexNeighbors,2)-innerSurf.vtcs(vertexList(iVertex),2);
  z = innerSurf.vtcs(vertexNeighbors,3)-innerSurf.vtcs(vertexList(iVertex),3);
  % get the rotation matrix that rotates coordinates so
  % that the normal at this vertex is [0 0 1];
  % first extract normal
  normal = vertexNormals(vertexList(iVertex),:);  
  normal = normal/norm(normal);
  % now projext normal on to YZ plane
  normalX = [normal(2:3)];normalX = normalX/norm(normalX);
  % calculate rotation to Z=1 vector in YZ plane
  thetaX = acos(dot(normalX,[0 1]));
  % and get the rotation matrix
  if ~isnan(thetaX)
    if normalX(1) < 0
      rotateX = [1 0 0;0 cos(pi-thetaX) sin(pi-thetaX);0 -sin(pi-thetaX) cos(pi-thetaX)]*diag([1 1 -1]);
    else
      rotateX = [1 0 0;0 cos(thetaX) sin(thetaX);0 -sin(thetaX) cos(thetaX)];
    end
    normal = normal*rotateX;
  else
    rotateX = eye(3);
  end
  % now do same for projection on to XZ plane
  normalY = [normal(1) normal(3)];normalY = normalY/norm(normalY);
  thetaY = acos(dot(normalY,[0 1]));
  if ~isnan(thetaY)
    if (normalY(1) < 0)
      rotateY = [cos(pi-thetaY) 0 sin(pi-thetaY);0 1 0;-sin(pi-thetaY) 0 cos(pi-thetaY)]*diag([1 1 -1]);
    else
      rotateY = [cos(thetaY) 0 sin(thetaY);0 1 0;-sin(thetaY) 0 cos(thetaY)];
    end     
  else
    rotateY = eye(3);
  end
  R = rotateX*rotateY;
  if 0
    normal = vertexNormals(vertexList(iVertex),:);  
    normal = normal/norm(normal);
    normal = normal*R;
    if (normal(1) > 0.1) || (normal(2) > 0.1) || (normal(3) < 0.9)
      disp(normal);
      keyboard
    end
  end
  % rotate the coordinates
  rotatedCoords = [x y z]*R;
  x = rotatedCoords(:,1);
  y = rotatedCoords(:,2);
  z = rotatedCoords(:,3);
  % now compute least-squares fit of a quadratic surface to the points
  quadraticCoef = pinv(0.5*[x.^2 2*x.*y y.^2])*z;
  % Now we make the cooeficients into a matrix
  A = [quadraticCoef(1) quadraticCoef(2);quadraticCoef(2) quadraticCoef(3)];
  % compute eigenvalues of this matrix. These are the principle curvatures;
  pCurvature = eig(A);
  % get the mean curvature
  m(iVertex) = mean(pCurvature);
  disppercent(iVertex/length(vertexList));
end
disppercent(inf);

% invert colors
m = -m;

% and smooth
for i=1:nSmooth
  m = connectionBasedSmooth(connectionMatrix,m);
end

if 0
surfname = fullfile(mrGetPref('volumeDirectory'),'jg19710606/jg041001ver2/jg_left_WM');
surf = loadSurfOFF(surfname);
surf.uniqueVertices = surf.vtcs;
surf.uniqueFaceIndexList = surf.tris;
connectionMatrix = double(findConnectionMatrix(surf));
d = dijkstra(connectionMatrix,10);

innerSurfname = fullfile(mrGetPref('volumeDirectory'),'jg19710606/jg041001ver2/jg_left_WM');
outerSurfname = fullfile(mrGetPref('volumeDirectory'),'jg19710606/jg041001ver2/jg_left_GM');
m = calcCurvature(innerSurfname,outerSurfname);
end

