function tseries = loadTSeries(view,scan,slice,frame)
%
% tSeries = loadTSeries(view,scan,[slice],[frame])
%
% Loads the tSeries corresponding to the specified scan and slice. The
% tseries files must be in <homedir>/<view>/<group>/TSeries, and must be
% nifti (or analyze) format (4D file containing all of the slices and
% temporal frames from a given scan).
%
% scan: number specfying which scan to load.
%
% slice: either a number specifying which slice to load or 'all' which
% loads all of the slices. Default: 'all'
%
% frame: optional number specifying which temporal frame to load. Ignored
% if slice is specified. Default: [] (load all frames).
%
% tseries: If a single slice is specified, tseries is returned as aa 2D
% matrix (nFrames x nVoxels). If slice is 'all' then the tseries is
% returned as a 4D array [x y z t]. If frame is specified, then the tseries
% is a single volume returned as a 3D array [x y t]
%
% djh 1/9/98
% djh 2/20/2001 Removed interpolation, dumped dat files
% JL 8/6/03 Updated to allow Analyze format.
% djh 5/2005 Updated to mrLoadRet 4.0

if ieNotDefined('slice')
	slice = 'all';
end
if ieNotDefined('frame')
    frame = [];
end

% Get the pathStr to the tseries file
pathStr = viewGet(view,'tseriesPathStr',scan);
% Error if file not found
if ~exist(pathStr,'file')
	mrErrorDlg(['File ',pathStr,' not found']);
end

% Load the tSeries
[path,name,ext,versn] = fileparts(pathStr);
if strcmp(slice,'all')
    if (isnumeric(frame) & (1 <= frame) & (frame <= viewGet(view,'nframes',scan)))
        % Load it
        [tseries,hdr] = cbiReadNifti(pathStr,{[],[],[],frame});
        dims = size(tseries);
    elseif isempty(frame)
        % Load it
        [tseries,hdr] = cbiReadNifti(pathStr,{[],[],[],frame});
        % Check # temporal frames
        dims = size(tseries);
        nFrames = dims(4);
        if (nFrames ~= viewGet(view,'totalFrames',scan))
            mrWarnDlg('loadTSeries: number of frames in tseries file does not match expected.');
        end
    else
        mrErrorDlg(['Invalid frame number: ',num2str(frame)]);
    end
    % Check data size
    dataSize = dims(1:3);
	if (dataSize ~= viewGet(view,'dataSize',scan))
		mrWarnDlg('loadTSeries: number of voxels in tseries does not match expected.');
	end
elseif (isnumeric(slice) & (1 <= slice) & (slice <= viewGet(view,'nslices',scan)))
	% Load single slice and reshape to nFrames X nVoxels
	[tseries,hdr] = cbiReadNifti(pathStr,{[],[],slice,[]});
	tseries = squeeze(tseries);
	% Reformat to nFrames x nVoxels
	dims = size(tseries);
	nFrames = dims(3);
	sliceDims = dims(1:2);
	nVoxels = prod(sliceDims);
	tseries = reshape(tseries,[nVoxels,nFrames])';
	% Check # voxels and # temporal frames
	if (nFrames ~= viewGet(view,'totalFrames',scan))
		mrWarnDlg('loadTSeries: number of frames in tseries file does not match expected.');
	end
	if sliceDims ~= viewGet(view,'sliceDims',scan)
		mrWarnDlg('loadTSeries: number of voxels in tseries does not match expected.');
	end
else
    mrErrorDlg(['Invalid slice number: ',slice]);
end
return
