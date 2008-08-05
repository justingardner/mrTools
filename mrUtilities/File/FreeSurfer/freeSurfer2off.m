% freeSurfer2off.m
%
%        $Id$
%      usage: freeSurfer2off(fsSurf, offSurf)
%         by: eli merriam
%       date: 07/11/07
%    purpose: 
%
function[] = freeSurfer2off(fsSurf, offSurf)

% check arguments
if (nargin ~= 2)
  help freeSurfer2off
  return
end

% read in the freesurfer file
[vertices, triangles] = freesurfer_read_surf(fsSurf);

% subtract 1 for OFF compatibility
triangles = triangles' -1;
vertices  = vertices'  -1;

% center image
vertices(1,:) = vertices(1,:) + 88;   % higher, more right
vertices(2,:) = vertices(2,:) + 128;    % higher, more anterior
vertices(3,:) = vertices(3,:) + 128;

% triangles(1) is number of vert/triangle: 3
% triangles(2:4) are the vertices of that triangles
% triangles(5) is color: 0
triangles = cat(1, repmat(3,1,length(triangles)), triangles, repmat(0,1,length(triangles)));

% write the OFF format file
fid = fopen(offSurf, 'w', 'ieee-be');

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [size(vertices,2) size(triangles,2) 0], 'int32'); 

% Vertices
fwrite(fid, vertices, 'float32');

% Faces
fwrite(fid, triangles, 'int32');

% Close file
fclose(fid);


