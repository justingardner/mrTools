function [view,filename] = saveNewTSeries(view,tseries,scanParams,hdr)
%
% view = saveNewTSeries(view,tseries,[scanParams],[hdr])
%
% tseries: tseries array (x,y,z,t) or a filename of a nifti file
%
% scanParams: structure with the following fields: fileName, description,
% junkFrames, nFrames. The other scanParams fields are set automatically
% from the nifti header. 
% Default:
%    scanParams.fileName = 'tseries-mmddyy-hhmmss.img'; 
%       where mmddyy = date and hhmmss = time
%    scanParams.junkFrames = 0;
%    scanParams.nFrames = size(tseries,4);
%    scanParams.description = '';
%
% Chooses between .img and .nii based on 'niftiFileExtension' preference.
%
% hdr: template for nifti header. The header is always passed through
% cbiCreateNiftiHeader to ensure consistency with the data. Default: [];
%
% $Id$
%

if ieNotDefined('scanParams')
  scanParams.fileName = [];
  scanParams.junkFrames = 0;
  scanParams.description = '';
end
if ieNotDefined('hdr')
  hdr = [];
end

if isempty(scanParams.fileName)
  scanParams.fileName = ['tseries-',datestr(now,'yymmdd-HHMMSS'),mrGetPref('niftiFileExtension')];
end
filename = scanParams.fileName;

% Save tseries 
tseriesdir = viewGet(view,'tseriesdir');
path = fullfile(tseriesdir,scanParams.fileName);

% see if the tseries is actually a string in which case we should
% copy the nifti file.
if isstr(tseries)
  success = copyNiftiFile(tseries,path);
  if ~success,return,end

  % get the number of frames
  hdr = cbiReadNiftiHeader(tseries);
  nFrames = hdr.dim(5);
else
  [byteswritten,hdr] = cbiWriteNifti(path,single(tseries),hdr,'float32');
  nFrames = size(tseries,4);
end

% set nFrames
if ~isfield(scanParams,'nFrames')
  scanParams.nFrames = nFrames;
end  
% Add new scan (note that the tseries file must already exist)
view = viewSet(view,'newScan',scanParams);
