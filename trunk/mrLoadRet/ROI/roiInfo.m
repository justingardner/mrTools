% roiInfo.m
%
%        $Id$
%      usage: roiInfo()
%         by: justin gardner
%       date: 01/16/08
%    purpose: 
%
function retval = roiInfo(v,roinum)

% check arguments
if ~any(nargin == [1 2])
  help roiInfo
  return
end

if ieNotDefined('roiNum')
  roiNum = viewGet(v,'currentROI');
end
roiName = viewGet(v,'roiName',roiNum);
roiDate = viewGet(v,'roidate',roiNum);
roiColor = viewGet(v,'roicolor',roiNum);
roiVoxelSize = viewGet(v,'roivoxelsize',roiNum);
roiVolume = viewGet(v,'roivolume',roiNum);
roiXform = viewGet(v,'roixform',roiNum);
roiNotes = viewGet(v,'roiNotes',roiNum);
vol2mag = viewGet(v,'roiVol2mag',roiNum);
vol2tal = viewGet(v,'roiVol2tal',roiNum);

% check to see which base anatomy this roi aligns with
baseMatch = {};
for bnum = 1:viewGet(v,'numberOfBaseVolumes')
  % get the base voxelSize and xfrom
  baseVoxelSize = viewGet(v,'baseVoxelSize',bnum);
  baseXform = viewGet(v,'baseXform',bnum);
  % if it matches, then put it in thee list of matching base names
  if isequal(baseXform,roiXform) && isequal(baseVoxelSize,roiVoxelSize)
    baseMatch{end+1} = viewGet(v,'baseName',bnum);
  end
end
if isempty(baseMatch),baseMatch = 'No matching base anatomy';,end
if length(baseMatch)==1,baseMatch = baseMatch{1};end

paramsInfo = {{'name',roiName,'editable=0','The name of the ROI'},...
  {'notes',roiNotes,'editable=0','Notes associated with ROI'},...
  {'date',roiDate,'editable=0','The date of creation'},...
  {'color',roiColor,'editable=0','ROI color'},...
  {'voxelsize',roiVoxelSize,'editable=0','Voxel dimensions in mm'},...
  {'volume',roiVolume,'editable=0','Volume of ROI in cubic mm'},...
  {'xform',roiXform,'editable=0','xform matrix specifies the transformation to the canoncial volume in magnet coordinates'},...
  {'vol2mag',vol2mag,'editable=0','xform matrix specifies the transformation of the canonical base volume to magnet coordinates'},...
  {'vol2tal',vol2tal,'editable=0','xform matrix specifies the transformation of the canonical base volume to talairach coordinates'},...
  {'baseMatch',baseMatch,'editable=0','The base volume that has the same voxel size and xform as this ROI. This is the base volume on which the ROI was originally defined. If there is no matching base anatomy, it means that the ROI was defined on a different base volume than the one you have loaded.'},...
  {'ROICoords',[],'type=pushbutton','buttonString=Show ROI coordinates','callback',@showCurrentROICoords,'callbackArg',v,'Print the coordinates for this ROI into the matlab window. Note that these will be the actual ROI coordinates not transformed into the scan coordinates. If you want the variable ROICoords set to the coordinates in your matlab workspace, you can hold the shift key down as you press this button (note that you have to have mgl in your path for this to work).'},...
  {'ROIScanCoords',[],'type=pushbutton','buttonString=Show scan coordinates','callback',@showCurrentROIScanCoords,'callbackArg',v,'Print the coordinates transformed into the scan coordinates for thie ROI to the matlab window. If you want the variable ROICoords set to the coordinates in your matlab workspace, you can hold the shift key down as you press this button (note that you have to have mgl in your path for this to work).'}};

% give ability to findROI for non baseCoordMapped ROIs
if isempty(viewGet(v,'baseCoordMap'))
  paramsInfo{end+1} = {'findROI',[],'type=pushbutton','buttonString=Find ROI','callback',@findROI,'callbackArg',v,'Go to the closest slice for which this ROI has some coordinates.'};
end

mrParamsDialog(paramsInfo,'ROI information',1.5);
