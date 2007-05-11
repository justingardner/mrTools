% NAME: mrAlign (version 4.3)
% AUTHOR: DJH
% DATE: 6/2004
% PURPOSE:
%	Aligning/registering inplane and volume anatomies.
% HISTORY:
%	Dozens of people contributed to the original mrAlign. It was written
%   using scripts in matlab 4 without any data structures. mrAlign4  
%   cleaned up the code. Version 4.1 adds a GUI. Version 4.2 ported to
%   matlab 7 and adds a number of additional features.


%%%%%%%%%%%%%%%%%%%
% Global Variables %
%%%%%%%%%%%%%%%%%%%

global ALIGN
global mrDEFAULTS

% Load .mrDefaults
mrDEFAULTS = loadMrDefaults;

% Check Matlab version number
[mlrVersion, expectedMatlabVersion] = mrLoadRetVersion;
version = ver('Matlab');
matlabVersion = str2num(version.Version(1:3));
if ~ismember(matlabVersion, expectedMatlabVersion);
    mrWarnDlg(['mrAlign is intended for Matlab ',num2str(expectedMatlabVersion),...
        '. You are running Matlab ',version.Version]);
end

ALIGN.version = mlrVersion;

ALIGN.volumePath = [];      % path string to volume anatomy file
ALIGN.volume = [];          % volume matrix
ALIGN.volumeHdr = [];       % loaded from nifti header
ALIGN.volSize = [];         % size of volume matrix
ALIGN.volumeVoxelSize = []; % 3x1 vector specifying voxel size
ALIGN.volumePermutation = eye(3);  % extracted from volume header to determine slice orientation
ALIGN.volumeClip = [];      % for display

ALIGN.inplanePath = [];     % path string to inplane anatomy file
ALIGN.inplanes = [];        % inplane matrix
ALIGN.inplaneHdr = [];      % loaded from nifti header
ALIGN.inplaneSize = [];     % size of inplane matrix
ALIGN.inplaneVoxelSize = [];% 3x1 vector specifying voxel size
ALIGN.inplanesClip = [];    % for display

ALIGN.crop = [];            % specifies 2D crop region used for contrast correction
ALIGN.guiXform = [];        % 4x4 transform matrix read from rot and trans GUI
ALIGN.xform = [];           % 4x4 transform matrix from inplane pixels -> volume pixels
ALIGN.NIter = 10;           % Number of iterations in auto alignment
ALIGN.sliceOrientation = 1; % 1=axial, 2=coronal, 3=sagittal
ALIGN.coords = [1,1,1];     % selected voxel

ALIGN.baseCmap = gray(256);
ALIGN.overlayCmap = hot(256);

disp(['mrAlign (version ',num2str(ALIGN.version),')']);

mrAlignGUI

return
