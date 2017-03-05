
%
%        $Id:$ 
%      usage: roi = mlrAnatDBGetPRFROI(v,roi)
%         by: justin gardner
%       date: 01/03/17
%    purpose: adds pRF parameters for voxels in ROI
%
%       e.g.: 
%
% roi = loadROITSeries(v,'lV1',1,'Concatenation','loadType=none');
% roi = mlrAnatDBGetPRFROI(v,roi);
%
function roi = mlrAnatDBGetPRFROI(v,roi,varargin)

% check arguments
if nargin < 2
  help mlrAnatDBGetPRFROI
  return
end

% parse arguments
getArgs(varargin,{'noPull=0'});

% check the roi
tf = isroi(roi);
if ~tf
  disp(sprintf('(mlrAnatDBGetPRFROI) Passed in argument is not a ROI'));
  return
end
if ~isfield(roi,'scanCoords')
  disp(sprintf('(mlrAnatDBGetPRFROI) roi does not have scanCoords - must be loaded from loadROITSeries'));
  return
end

% get the pRF for this subjet id
subjectID = roi.subjectID;
if isempty(subjectID),subjectID = v;end
pRF = mlrAnatDBGetPRF(subjectID,'noPull',noPull);
if isempty(pRF)
  disp(sprintf('(mlrAnatDBGetPRFROI) Could not load pRF'));
  return
end

% check xforms
if isempty(roi.vol2mag) || isempty(roi.scan2roi)
  disp(sprintf('(mlrAnatDBGetPRFROI) Empty vol2mag or scan2mag on roi, cannot xform coordinates'));
  return
end

% get relevant d from pRF
d = pRF.d{pRF.scanNum};

% transform roi coordinates into pRF coordinates
scanCoords = roi.scanCoords;
scanCoords(4,:) = 1;
roipRFCoords = round(inv(pRF.scan2mag)*roi.vol2mag*roi.scan2roi*scanCoords);
roipRFLinearCoords = mrSub2ind(pRF.scanDims,roipRFCoords(1,:),roipRFCoords(2,:),roipRFCoords(3,:));
[dummy roiMatch pRFMatch] = intersect(roipRFLinearCoords,d.linearCoords);

% get the params info struct for interpreting parameters
roi.pRF.paramsInfo = d.paramsInfo;

% get parameters
nParams = length(d.paramsInfo.paramNames);
roi.pRF.params = nan(nParams,roi.n);
roi.pRF.params(:,roiMatch) = d.params(:,pRFMatch);

% get r
roi.pRF.r = nan(roi.n,1);
roi.pRF.r(roiMatch) = d.r(pRFMatch);

% debugging code - display pRF
if 0
mlrSmartfig('mlrAnatDBGetPRFROI','reuse');clf;
hold on
for i = 1:roi.n
  % get pRF coords
  x = roi.pRF.params(1,i);
  y = roi.pRF.params(2,i);
  sigma = roi.pRF.params(3,i);
  if ~isnan(x) && ~isnan(y) && ~isnan(sigma) && (sigma>0)
    h = rectangle('Position',[x y sigma sigma],'Curvature',[1 1]);
  end
end
end
