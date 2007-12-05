% writeOFF.m
%
%       $Id$	
%      usage: writeOFF(surf, outName)
%         by: eli merriam
%       date: 10/25/07
%    purpose: 
%
function retval = writeOFF(surf, outName)

% write the OFF format file

% undo what loadSurfOFF did in loading surface into matlab
surf.tris = surf.tris - 1;

surf.tris = cat(1, repmat(3,1,length(surf.tris)), surf.tris', repmat(0,1,length(surf.tris)));

surf.vtcs = [surf.vtcs(:,2) surf.vtcs(:,1) surf.vtcs(:,3)];
surf.vtcs = surf.vtcs' - 2;


fid = fopen(outName, 'w', 'ieee-be');

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [surf.Nvtcs surf.Ntris 0], 'int32'); 

% Vertices
fwrite(fid, surf.vtcs, 'float32');

% Faces
fwrite(fid, surf.tris, 'int32');

% Close file
fclose(fid);

