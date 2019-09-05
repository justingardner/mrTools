% mlrTestROIsInMLR.m
%
%      usage: mlrTestROIsInMLR(rois,surfaces)
%         by: justin gardner
%       date: 09/05/19
%    purpose: Creates a temporary MLR directory and show a set of ROIs in there
%             for use with mlrImportNeuropythy and mlrImportFreesurferLabels
%
function retval = mlrTestROIsInMLR(rois,surfaces)

% check arguments
if nargin < 2
  help mlrTestROIsInMLR
  return
end

% try out in an empty MLR directory
mrQuit(0);
testDirname = sprintf('mlrTestROIs_%s%s',datestr(now,'YYMMDD'),datestr(now,'HHMMSS'));
makeEmptyMLRDir(testDirname,'defaultParams=1');
curpwd = pwd;
cd(testDirname)

% open up the view
mrLoadRet;

% load the surfaces
for iSurface = 1:length(surfaces)
  base = importSurfaceOFF(surfaces{iSurface});
  viewSet(getMLRView,'newBase',base);
end

% load the rois
for iROI = 1:length(rois)
  viewSet(getMLRView,'newROI',rois{iROI});
  viewSet(getMLRView,'currentROI',viewGet(getMLRView,'nrois'));
end

% refresh the display
refreshMLRDisplay(viewGet(getMLRView,'viewNum'));

% wait till user does a dbcont
dispHeader(sprintf('(mlrTestROIsInMLR) Type dbcont when you are done checking ROIs in MLR'));
keyboard

% quit and remove the temporary MLR directory
mrQuit(0);
cd(curpwd);
system(sprintf('rm -rf %s',testDirname));

