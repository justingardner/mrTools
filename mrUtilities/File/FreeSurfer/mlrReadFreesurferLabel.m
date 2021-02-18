function [vertices,header] = mlrReadFreesurferLabel(filename)

% open file
fileID = fopen(filename);
if fileID == -1
  mrErrorDlg(sprintf('(mlrReadFreesurferLabel) Could not open file %s',filename));
end

% get header
header = fgetl(fileID);

% get number of vertices
nVertices = str2num(fgetl(fileID));
% get each vertex row-by-row
for iVertex = 1:nVertices
  nextLine = fgetl(fileID);
  if ischar(nextLine)
    vertices(iVertex,:) = str2num(nextLine);
  else
    fclose(fileID);
    mrErrorDlg(sprintf('(mlrReadFreesurferLabel) There was an issue reading line %d of label file %s',iVertex,filename));
  end
end

% close file
fclose(fileID);