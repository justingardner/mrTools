function eventRelatedPlot(view,overlayNum,scan,x,y,s)
% eventRelatedPlot.m
%
%      usage: eventRelatedPlot()
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%


% check arguments
if ~any(nargin == [1:6])
  help eventRelatedPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};

% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot');

% get the number of ROIs, we first will 
% look to see if the user clicked on a 
% voxel that is contained in an ROI.
roi = {};
nROIs = viewGet(view,'numberOfROIs');
for roinum = 1:nROIs
  roicoords = getRoiCoordinates(view,scan,roinum);
  % see if this is a matching roi
  if ismember([x y s],roicoords','rows')
    % get the roi
    roi{end+1} = viewGet(view,'roi',roinum);
    % change the coordinates to our coordinates
    roi{end}.coords = roicoords;
  end
end

if isempty(d)
  disp('No analysis');
  reutrn
end
 
% get the estimated hemodynamic responses
[ehdr time] = gethdr(d,x,y,s);

% plot the timecourse
subplot(2,2,1:2)
tSeries = squeeze(loadTSeries(view,scan,s,[],x,y));
plot(tSeries);
xlabel('Volume number');
ylabel('MRI signal');
% and the stimulus times
hold on
for i = 1:d.nhdr
  vline(d.stimvol{i},getcolor(i));
end
axis tight;
subplot(2,2,3);
% and display ehdr
plotEhdr(time,ehdr);
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',x,y,s,analysis.overlays.data{scan}(x,y,s)));
xaxis(0,d.hdrlen*d.tr);
if isfield(d,'stimNames')
  legend(d.stimNames);
end

% if there is a roi, compute its average hemodynamic response
for roinum = 1:length(roi)
  subplot(2,2,4);
  ehdr = [];
  for voxnum = 1:size(roi{roinum}.coords,2)
    [ehdr(voxnum,:,:) time] = gethdr(d,roi{roinum}.coords(1:3,voxnum));
  end
  plotEhdr(time,squeeze(mean(ehdr)));
  title(sprintf('%s (n=%i)',roi{roinum}.name,size(roi{roinum}.coords,2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr)

% and display ehdr
for i = 1:size(ehdr,1)
  h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,'-')),'MarkerSize',8);
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the roi coordinates in this overlay coordinate frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overlayCoords = getRoiCoordinates(view,scan,roinum)

% get the base transform
baseNum = viewGet(view,'currentBase');
baseXform = viewGet(view,'baseXform',baseNum);
baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);

% get the base to overlay scale factors
baseDims = viewGet(view,'baseDims',baseNum);
overlayDims = viewGet(view,'overlayDims',scan);

% Get subset of coords corresponding to the current slice
roiCoords = viewGet(view,'roiCoords',roinum);
roiXform = viewGet(view,'roiXform',roinum);
roiVoxelSize = viewGet(view,'roiVoxelSize',roinum);

% Use xformROI to supersample the coordinates
baseCoords = round(xformROIcoords(roiCoords,inv(baseXform)*roiXform,roiVoxelSize,baseVoxelSize));

% and convert them into the overlay units by noting dimension
% scale factor (is this the best way to do this?)
if ~isempty(baseCoords)
  overlayCoords(1,:) = round(baseCoords(1,:)*overlayDims(1)/baseDims(1));
  overlayCoords(2,:) = round(baseCoords(2,:)*overlayDims(2)/baseDims(2));
  overlayCoords(3,:) = baseCoords(3,:);
else
  overlayCoords = [];
end
% return the unique ones
overlayCoords = unique(overlayCoords','rows')';