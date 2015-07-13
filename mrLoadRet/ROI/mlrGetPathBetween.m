% mlrGetPathBetween.m
%
%        $Id:$ 
%      usage: pathList = mlrGetPathBetween(tris,v1,v2)
%         by: justin gardner
%       date: 07/12/15
%    purpose: Gets shortest path between two vertices given triangle list
%
%       e.g.: baseSurface = viewGet(v,'baseSurface');
%             mlrGetPathBetween(baseSurface.tris, 100,150);
%
function pathList = mlrGetPathBetween(tris,v1,v2)

% check arguments
if ~any(nargin == [3])
  help mlrGetPathBetween
  return
end

% start with v1 in edge list
edgeList = v1;

% start the predecessor list with v1
predList = v1;

% get surface
trisSize = size(tris);

% keep going until we have v1 in our edgeList
while ~ismember(v2,edgeList)
  % find all the triangles that have one vertex in common with the current edgeList of vertices
  [intersectTris predecessors] = ismember(tris(:),edgeList);
  intersectTris = find(intersectTris);
  % keep also what vertex from the edgeList each of these comes from
  predecessors = edgeList(predecessors(find(predecessors)));
  % now find all the triangles that they intersect with - this means find 
  % the triangle number (since the above might have made a match at any one of
  % the three vertices in a triangle).
  [intersectTris,dummy] = ind2sub(trisSize,intersectTris);
  [intersectTris uniqueIndexes] = unique(intersectTris);
  predecessors = predecessors(uniqueIndexes);
  % get al vertices from the triangles found before (this is now the list of new vertices)
  [newEdgeList uniqueIndex] = unique(tris(sub2ind(trisSize,repmat(intersectTris,1,3),repmat([1 2 3],length(intersectTris),1))));
  predecessors = repmat(predecessors,1,3);
  predecessors = predecessors(uniqueIndex);
  % find all the vertices that were in the original
  [oldEdges oldEdgeIndex] = ismember(newEdgeList,edgeList);
  % and assign those to have the predecessor list from before
  predecessors(oldEdges) = predList(oldEdgeIndex(oldEdges));
  predList = predecessors;
  edgeList = newEdgeList;
end

% now go backwords from v2 and find the path to v1
pathList = v2;
while pathList(end) ~= v1
  pathList(end+1) = predList(find(edgeList == pathList(end)));
end


