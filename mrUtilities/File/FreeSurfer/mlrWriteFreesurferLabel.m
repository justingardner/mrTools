function mlrWriteFreesurferLabel(filename,header,vertices)

% open file
 fileID = fopen(filename,'w');
if fileID == -1
  mrErrorDlg(sprintf('(mlrReadFreesurferLabel) Could not open file %s for writing',filename));
end

fprintf(fileID,header);
fprintf(fileID,'%d\n', size(vertices,1));
fprintf(fileID,'%d  %.3f  %.3f  %.3f %.10f\n', vertices');
fclose(fileID);
  
