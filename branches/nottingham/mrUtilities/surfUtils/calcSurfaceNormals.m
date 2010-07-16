% calcNormals.m
%
%        $Id$ 
%      usage: normals = calcNormals(surf)
%         by: justin gardner
%       date: 08/10/08
%    purpose: Calculates vertex normals from a surf returned by loadSurfOFF
%             This is computes by averaging the surface normals of each triangle
%             a vertex is a part of
function vertexNormals = calcSurfaceNormals(surf)

% check arguments
if ~any(nargin == [1])
  help calcNormals
  return
end

% first compute the normals to each triangle face.
% this is done with the cross product of two edge vectors
triNormals = zeros(surf.Ntris,3);
disppercent(-inf,'(calcSurfaceNormal) Computing triangle normals');
for iTri = 1:surf.Ntris
  % get the three vertices of this triangle
  vertex1 = surf.vtcs(surf.tris(iTri,1),:);
  vertex2 = surf.vtcs(surf.tris(iTri,2),:);
  vertex3 = surf.vtcs(surf.tris(iTri,3),:);
  % and compute the surface normal using the cross product
  triNormals(iTri,:) = cross(vertex2-vertex1,vertex2-vertex3);
  triNormals(iTri,:) = triNormals(iTri,:)/norm(triNormals(iTri,:));
  disppercent(iTri/surf.Ntris);
end
disppercent(inf);

% now compute the normal to each vertex as the average
% normal of all the triangle faces it belongs to
disppercent(-inf,'(calcSurfaceNormal) Computing vertex normals');
vertexNormals = zeros(surf.Nvtcs,3);
for iVtx = 1:surf.Nvtcs
  % get which triangles this vertex belongs to
  [triNums edgeNums] = find(iVtx == surf.tris);
  % and then get the mean of the normals to those triangls
  vertexNormals(iVtx,:) = mean(triNormals(triNums,:));
  disppercent(iVtx/surf.Nvtcs);
end
disppercent(inf);
