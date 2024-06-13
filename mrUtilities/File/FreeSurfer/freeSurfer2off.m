% freeSurfer2off.m
%
%        $Id$
%      usage: freeSurfer2off(fsSurf, offSurf, <volumeSize>, <volumeSize>)
%         by: eli merriam
%       date: 07/11/07
%    purpose: Converts vertices from free surfer conventions and saves as an off. Note that
%             freeSurfer is 1 based and coordinates start in the middle of the volume. We
%             therefore have to add half the volume size (in mm) to the coordinates to convert.
%             The default is to assume that the volumeSize is 176x256x256 and the pixelSize 1x1x1.
%             Note that the script mlrImportFreeSurfer crops volumes to that size.
%
function [] = freeSurfer2off(fsSurf, offSurf, volumeSize, pixelSize)

% check arguments
if (nargin < 2)
  help freeSurfer2off
  return
end

% default volume size in RAS coordinates
if nargin < 3
  volumeSize = [176 256 256];
end
if nargin < 4
  pixelSize = [1 1 1];
end
  

% read in the freesurfer file
[vertices, triangles] = freesurfer_read_surf(fsSurf);

% subtract 1 for OFF compatibility
triangles = triangles' -1; % JB (01/08/2020): this -1 is related to 0-indexing used in the surfRelax format.
vertices  = vertices'  -1; % I dont' think this -1 represents the same thing because it is subtracted from coordinates in mm
                           % Instead, I think it's necessary because Freesurfer (mri_convert) always puts ones in hdr.pixdim(6:8) and
                           % surfRelax assumes (maybe wrongly?) that these represent x,y,z offsets and uses them to convert the coordinates to "world2 coordinates"
                           % This is then undone when converting back to array coordinates using xformSurfaceWorld2Array or mlrXFormFromHeader(hdr,'world2array')
                           % (but this time using the actual values in hdr.pixdim(6:8), which are always 1
% center image
vertices(1,:) = vertices(1,:) + pixelSize(1)*volumeSize(1)/2;   % higher, more right
vertices(2,:) = vertices(2,:) + pixelSize(2)*volumeSize(2)/2;   % higher, more anterior
vertices(3,:) = vertices(3,:) + pixelSize(3)*volumeSize(3)/2;   % higher, more superior

% triangles(1) is number of vert/triangle: 3
% triangles(2:4) are the vertices of that triangles
% triangles(5) is color: 0
triangles = cat(1, repmat(3,1,length(triangles)), triangles, zeros(1,length(triangles)));

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


