function [view,filename] = saveNewTSeries(view,tseries,scanParams,hdr)
%
% view = saveNewTSeries(view,tseries,[scanParams],[hdr])
%
% tseries: tseries array (x,y,z,t)
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
	scanParams.nFrames = size(tseries,4);
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
[byteswritten,hdr] = cbiWriteNifti(path,tseries,hdr);

% Add new scan (note that the tseries file must already exist)
view = viewSet(view,'newScan',scanParams);
