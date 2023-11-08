%
% function thisView = makeSphereROIatCoords(thisView,params,<'justGetParams'>)
%
%   Purpose: make a spherical ROI centered at specified coordinates and add it to the mrLoadRet view
%   
%   Usage:   [~,params] = makeSphereROIatCoords(thisView,[],'justGetParams') % to get default parameters
%            % then set the parameters and use them to create the ROI in the desired location, e.g.:
%            params.centerCoords = [-45 -57 -12]; % specify the coordinates in MNI coordinates (default)
%            params.name = 'myROI'; % specify the ROI name
%            thisView = makeSphereROIatCoords(thisView,params); % make the sphere ROI and add it to the view
%            refreshMLRDisplay(thisView); % display the ROI in the GUI
%
%   Params: - name: name of the ROI. if left empty, the name will be ROI1, ROI2, etc..., whichever number
%                   if available (not yet used in the view)
%           - color: color of the ROI
%           - centerCoords: coordinates of the center of the sphere (see coordSpace for the coordinate system used)
%           - coordsSpace: space in which the coordinates are specified: options are 'base': the current base
%                          or 'MNI' (default): MNI space (mniInfo must be defined, see mlrImportSPMnormalization)
%           - radius: radius of the ROI in mm (in magnet space, i.e. distances in the participant's space)
%

function [thisView, params] = makeSphereROIatCoords(thisView,params,varargin)

eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
% if ieNotDefined('defaultParams'),defaultParams = 0;end

if ~ismember(nargin,[1 2 3])
  help('makeSphereROIatCoords');
  return
end

if ieNotDefined('thisView')
  mrWarnDlg('(makeSphereROIatCoords) No valid mrLoadRet view was provided');
  return
end

% set default parameters
if ieNotDefined('params')
  params = struct;
end

if fieldIsNotDefined(params,'name')
  params.name = ''; % name of the ROI. if left empty, the name will be ROI1, ROI2, etc..., whichever number if available (not yet used in the view)
end
if fieldIsNotDefined(params,'color')
  params.color = 'black'; % color of the ROI
end
if fieldIsNotDefined(params,'centerCoords')
  params.centerCoords = [0 0 0]; % coordinates of the center of the sphere (see coordSpace for the coordinate system used)
end
if fieldIsNotDefined(params,'coordsSpace')
  params.coordsSpace = 'MNI'; % space in which the coordinates are specified: options are 'base': the current base (default) or 'MNI': MNI space (mniInfo must be defined, see mlrImportSPMnormalization)
end
if fieldIsNotDefined(params,'radius')
  params.radius = 10; % radius of the ROI in mm (in magnet space, i.e. distances in the participant's space)
end

if justGetParams
  return;
end

base2mag = viewGet(thisView,'base2mag');
switch lower(params.coordsSpace)
  case 'mni'
    mniInfo = viewGet(thisView,'mniInfo');
    if isempty(mniInfo)
      mrWarnDlg('(makeSphereROIatCoords) MNI normalization info is not defined for this participant/session. You must first run mlrImportSPMnormalization.')
      return;
    end
    
    mniCoords = params.centerCoords;
    mniVolCoords = mniInfo.mni2mnivol*[mniCoords 1]';
    
    [mniVolGridX,mniVolGridY,mniVolGridZ] = ndgrid(1:size(mniInfo.mnivol2magCoordMap,1),1:size(mniInfo.mnivol2magCoordMap,2),1:size(mniInfo.mnivol2magCoordMap,3));
    magCenterCoords(1) = interpn(mniVolGridX,mniVolGridY,mniVolGridZ, mniInfo.mnivol2magCoordMap(:,:,:,1), mniVolCoords(1), mniVolCoords(2), mniVolCoords(3));
    magCenterCoords(2) = interpn(mniVolGridX,mniVolGridY,mniVolGridZ, mniInfo.mnivol2magCoordMap(:,:,:,2), mniVolCoords(1), mniVolCoords(2), mniVolCoords(3));
    magCenterCoords(3) = interpn(mniVolGridX,mniVolGridY,mniVolGridZ, mniInfo.mnivol2magCoordMap(:,:,:,3), mniVolCoords(1), mniVolCoords(2), mniVolCoords(3));
    %   coords.mag.x = mniInfo.mnivol2magCoordMap(round(mniVolCoords(1)),round(mniVolCoords(2)),round(mniVolCoords(3)),1); %
    %   coords.mag.y = mniInfo.mnivol2magCoordMap(round(mniVolCoords(1)),round(mniVolCoords(2)),round(mniVolCoords(3)),2); % less accurate but faster?
    %   coords.mag.z = mniInfo.mnivol2magCoordMap(round(mniVolCoords(1)),round(mniVolCoords(2)),round(mniVolCoords(3)),3); %
    
  case 'base'
    magCenterCoords = base2mag * [params.centerCoords 1]';

end

roi = makeEmptyROI(thisView, 'scanNum', 0, 'name',params.name);
roi.color = params.color;

% find all coordinates within radius of center coordinates
switch(viewGet(thisView,'baseType'))
  case 0
    baseDims = viewGet(thisView,'basedims');
  case {1,2}
    baseCoordMap = viewGet(thisView,'baseCoordMap');
    baseDims = baseCoordMap.dims;
end

[allBaseCoordsX,allBaseCoordsY,allBaseCoordsZ] = ndgrid(1:baseDims(1),1:baseDims(2),1:baseDims(3));
allMagCoords = base2mag * [allBaseCoordsX(:) allBaseCoordsY(:) allBaseCoordsZ(:) ones(prod(baseDims),1)]';
whichCoords = sqrt((allMagCoords(1,:) - magCenterCoords(1)).^2 + (allMagCoords(2,:)- magCenterCoords(2)).^2 +(allMagCoords(3,:) - magCenterCoords(3)).^2 ) < params.radius;
roi.coords = unique(round(base2mag \ allMagCoords(:,whichCoords))','rows')';

if isempty(roi.coords)
  mrWarnDlg('(makeSphereROIatCoords) The sphere has no coordinates in the current base volume. No ROI was added to the view.');
  return;
end

% Add it to the view
thisView= viewSet(thisView,'newROI',roi);

