% writeOFF.m
%
%      usage: writeOFF()
%         by: eli merriam
%       date: 10/25/07
%    purpose: 
%
function retval = writeOFF(surf, outName)
% write the OFF format file

fid = fopen(outName, 'w', 'ieee-be');

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [surf.Nvtcs surf.Ntris 0], 'int32'); 

% Vertices
fwrite(fid, surf.vtcs, 'float32');

% Faces
fwrite(fid, surf.tris, 'int32');

% Close file
fclose(fid);

